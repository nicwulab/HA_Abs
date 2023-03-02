import tensorflow as tf
import numpy as np

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Concatenate, Embedding, Dense, Dropout, Input, Layer, LayerNormalization, GlobalAveragePooling1D
from tensorflow.keras import Model
# from tensorflow.keras.utils import to_categorical
# from tensorflow.keras import layers

def get_angles(pos, i, d_model):
    angle_rates = 1 / np.power(10000, (2 * (i//2)) / np.float32(d_model))
    return pos * angle_rates

def positional_encoding(position, d_model):
    angle_rads = get_angles(np.arange(position)[:, np.newaxis],
                            np.arange(d_model)[np.newaxis, :],
                            d_model)

    # apply sin to even indices in the array; 2i
    angle_rads[:, 0::2] = np.sin(angle_rads[:, 0::2])

    # apply cos to odd indices in the array; 2i+1
    angle_rads[:, 1::2] = np.cos(angle_rads[:, 1::2])

    pos_encoding = angle_rads[np.newaxis, ...]

    return tf.cast(pos_encoding, dtype=tf.float32)

def create_padding_mask(seq):
    seq = tf.cast(tf.math.equal(seq, 0), tf.float32)

    # add extra dimensions to add the padding
    # to the attention logits.
    return seq[:, tf.newaxis, tf.newaxis, :]  # (batch_size, 1, 1, seq_len)

def create_look_ahead_mask(size):
    mask = 1 - tf.linalg.band_part(tf.ones((size, size)), -1, 0)
    return mask  # (seq_len, seq_len)

def scaled_dot_product_attention(q, k, v, mask):
    """Calculate the attention weights.
    q, k, v must have matching leading dimensions.
    k, v must have matching penultimate dimension, i.e.: seq_len_k = seq_len_v.
    The mask has different shapes depending on its type(padding or look ahead)
    but it must be broadcastable for addition.

    Args:
      q: query shape == (..., seq_len_q, depth)
      k: key shape == (..., seq_len_k, depth)
      v: value shape == (..., seq_len_v, depth_v)
      mask: Float tensor with shape broadcastable
            to (..., seq_len_q, seq_len_k). Defaults to None.

    Returns:
      output, attention_weights
    """

    # (..., seq_len_q, seq_len_k)
    matmul_qk = tf.matmul(q, k, transpose_b=True)

    # scale matmul_qk
    dk = tf.cast(tf.shape(k)[-1], tf.float32)
    scaled_attention_logits = matmul_qk / tf.math.sqrt(dk)

    # add the mask to the scaled tensor.
    if mask is not None:
        scaled_attention_logits += (mask * -1e9)

    # softmax is normalized on the last axis (seq_len_k) so that the scores
    # add up to 1.
    attention_weights = tf.nn.softmax(
        scaled_attention_logits, axis=-1)  # (..., seq_len_q, seq_len_k)

    output = tf.matmul(attention_weights, v)  # (..., seq_len_q, depth_v)

    return output, attention_weights

def print_out(q, k, v):
    temp_out, temp_attn = scaled_dot_product_attention(
        q, k, v, None)
    print('Attention weights are:')
    print(temp_attn)
    print('Output is:')
    print(temp_out)

class Defined_MultiHeadAttention(Layer):
    def __init__(self, d_model, num_heads):
        super(Defined_MultiHeadAttention, self).__init__()
        self.num_heads = num_heads
        self.d_model = d_model

        assert d_model % self.num_heads == 0

        self.depth = d_model // self.num_heads

        self.wq = Dense(d_model)
        self.wk = Dense(d_model)
        self.wv = Dense(d_model)

        self.dense = Dense(d_model)

    def split_heads(self, x, batch_size):
        """Split the last dimension into (num_heads, depth).
        Transpose the result such that the shape is (batch_size, num_heads, seq_len, depth)
        """
        x = tf.reshape(x, (batch_size, -1, self.num_heads, self.depth))
        return tf.transpose(x, perm=[0, 2, 1, 3])

    def call(self, v, k, q, mask):
        batch_size = tf.shape(q)[0]

        q = self.wq(q)  # (batch_size, seq_len, d_model)
        k = self.wk(k)  # (batch_size, seq_len, d_model)
        v = self.wv(v)  # (batch_size, seq_len, d_model)

        # (batch_size, num_heads, seq_len_q, depth)
        q = self.split_heads(q, batch_size)
        # (batch_size, num_heads, seq_len_k, depth)
        k = self.split_heads(k, batch_size)
        # (batch_size, num_heads, seq_len_v, depth)
        v = self.split_heads(v, batch_size)

        # scaled_attention.shape == (batch_size, num_heads, seq_len_q, depth)
        # attention_weights.shape == (batch_size, num_heads, seq_len_q, seq_len_k)
        scaled_attention, attention_weights = scaled_dot_product_attention(
            q, k, v, mask)

        # (batch_size, seq_len_q, num_heads, depth)
        scaled_attention = tf.transpose(scaled_attention, perm=[0, 2, 1, 3])

        concat_attention = tf.reshape(scaled_attention,
                                      (batch_size, -1, self.d_model))  # (batch_size, seq_len_q, d_model)

        # (batch_size, seq_len_q, d_model)
        output = self.dense(concat_attention)

        return output, attention_weights

def point_wise_feed_forward_network(d_model, dff):
    # return tf.keras.Sequential([
    return Sequential([
        # (batch_size, seq_len, dff)
        Dense(dff, activation='relu'),
        Dense(d_model)  # (batch_size, seq_len, d_model)
    ])

class EncoderLayer(Layer):
    def __init__(self, d_model, num_heads, dff, rate=0.1):
        super(EncoderLayer, self).__init__()

        self.mha = Defined_MultiHeadAttention(d_model, num_heads)
        self.ffn = point_wise_feed_forward_network(d_model, dff)

        self.layernorm1 = LayerNormalization(epsilon=1e-6)
        self.layernorm2 = LayerNormalization(epsilon=1e-6)

        self.dropout1 = Dropout(rate)
        self.dropout2 = Dropout(rate)

    def call(self, x, training, mask):

        # (batch_size, input_seq_len, d_model)
        attn_output, _ = self.mha(x, x, x, mask)
        attn_output = self.dropout1(attn_output, training=training)
        # (batch_size, input_seq_len, d_model)
        out1 = self.layernorm1(x + attn_output)

        ffn_output = self.ffn(out1)  # (batch_size, input_seq_len, d_model)
        ffn_output = self.dropout2(ffn_output, training=training)
        # (batch_size, input_seq_len, d_model)
        out2 = self.layernorm2(out1 + ffn_output)

        return out2

class DecoderLayer(Layer):
    def __init__(self, d_model, num_heads, dff, rate=0.1):
        super(DecoderLayer, self).__init__()

        self.mha1 = Defined_MultiHeadAttention(d_model, num_heads)
        self.mha2 = Defined_MultiHeadAttention(d_model, num_heads)

        self.ffn = point_wise_feed_forward_network(d_model, dff)

        self.layernorm1 = LayerNormalization(epsilon=1e-6)
        self.layernorm2 = LayerNormalization(epsilon=1e-6)
        self.layernorm3 = LayerNormalization(epsilon=1e-6)

        self.dropout1 = Dropout(rate)
        self.dropout2 = Dropout(rate)
        self.dropout3 = Dropout(rate)

    def call(self, x, enc_output, training,
             look_ahead_mask, padding_mask):
        # enc_output.shape == (batch_size, input_seq_len, d_model)

        # (batch_size, target_seq_len, d_model)
        attn1, attn_weights_block1 = self.mha1(x, x, x, look_ahead_mask)
        attn1 = self.dropout1(attn1, training=training)
        out1 = self.layernorm1(attn1 + x)

        attn2, attn_weights_block2 = self.mha2(
            enc_output, enc_output, out1, padding_mask)  # (batch_size, target_seq_len, d_model)
        attn2 = self.dropout2(attn2, training=training)
        # (batch_size, target_seq_len, d_model)
        out2 = self.layernorm2(attn2 + out1)

        ffn_output = self.ffn(out2)  # (batch_size, target_seq_len, d_model)
        ffn_output = self.dropout3(ffn_output, training=training)
        # (batch_size, target_seq_len, d_model)
        out3 = self.layernorm3(ffn_output + out2)

        return out3, attn_weights_block1, attn_weights_block2

class Encoder(Layer):
    def __init__(self, num_layers, d_model, num_heads, dff, input_vocab_size,
                 maximum_position_encoding, rate=0.1):
        super(Encoder, self).__init__()

        self.d_model = d_model
        self.num_layers = num_layers

        self.embedding = Embedding(input_vocab_size, d_model)
        self.pos_encoding = positional_encoding(maximum_position_encoding,
                                                self.d_model)

        self.enc_layers = [EncoderLayer(d_model, num_heads, dff, rate)
                           for _ in range(num_layers)]

        self.dropout = Dropout(rate)

    def call(self, x, training, mask):

        seq_len = tf.shape(x)[1]

        # adding embedding and position encoding.
        x = self.embedding(x)  # (batch_size, input_seq_len, d_model)
        x *= tf.math.sqrt(tf.cast(self.d_model, tf.float32))
        x += self.pos_encoding[:, :seq_len, :]

        x = self.dropout(x, training=training)

        for i in range(self.num_layers):
            x = self.enc_layers[i](x, training, mask)

        return x  # (batch_size, input_seq_len, d_model)

def create_mlp(dim=30):
    # define our MLP network
    model = Sequential()
    model.add(Dense(dim, activation="relu"))
    model.add(Dense(dim, activation="relu"))

    return model

def build_t_encoder(
    num_layers,
    d_model,
    num_heads,
    dff,
    input_vocab_size,
    maximum_position_encoding,
    training=True,
    rate=0.1,
    shape=(30,)
):
    inputs = Input(shape=shape)
    x=inputs
    enc_padding_mask = create_padding_mask(inputs)
    encoder = Encoder(num_layers, d_model, num_heads, dff,input_vocab_size, maximum_position_encoding, rate)
    x = encoder(x,training,enc_padding_mask)
    x = GlobalAveragePooling1D(data_format="channels_first")(x)
    t_encoder = Model(inputs, x)
    # for dim in mlp_units:
    #     x = Dense(dim, activation="relu")(x)
    #     x = layers.Dropout(rate)(x)
    # outputs = Dense(n_classes, activation="softmax")(x)
    return t_encoder

# from options import opt
# multi_class = opt.multi_class

def CDR_model_original6(max_length, mlp_units=[512,128,64],n_classes=1,rate=0.1, multi_class=False):
    shapes = [(m,) for m in max_length]
    input1=Input(shape=shapes[0])
    input2=Input(shape=shapes[1])
    input3=Input(shape=shapes[2])
    input4=Input(shape=shapes[3])
    input5=Input(shape=shapes[4])
    input6=Input(shape=shapes[5])

    # [print(x.shape) for x in [input1, input2, input3, input4, input5, input6]]

    encoder1=build_t_encoder(num_layers=4,d_model=256,num_heads=4,dff=512,input_vocab_size=22,maximum_position_encoding=max_length[0],shape=shapes[0])(input1)
    encoder2=build_t_encoder(num_layers=4,d_model=256,num_heads=4,dff=512,input_vocab_size=22,maximum_position_encoding=max_length[1],shape=shapes[1])(input2)
    encoder3=build_t_encoder(num_layers=4,d_model=256,num_heads=4,dff=512,input_vocab_size=22,maximum_position_encoding=max_length[2],shape=shapes[2])(input3)
    encoder4=build_t_encoder(num_layers=4,d_model=256,num_heads=4,dff=512,input_vocab_size=22,maximum_position_encoding=max_length[3],shape=shapes[3])(input4)
    encoder5=build_t_encoder(num_layers=4,d_model=256,num_heads=4,dff=512,input_vocab_size=22,maximum_position_encoding=max_length[4],shape=shapes[4])(input5)
    encoder6=build_t_encoder(num_layers=4,d_model=256,num_heads=4,dff=512,input_vocab_size=22,maximum_position_encoding=max_length[5],shape=shapes[5])(input6)

    x = Concatenate()([encoder1, encoder2,encoder3, encoder4,encoder5, encoder6])
    for dim in mlp_units:
        x = Dense(dim, activation="relu")(x)
        x = Dropout(rate)(x)
    act = 'softmax' if multi_class else 'sigmoid'
    outputs = Dense(n_classes, activation=act)(x)
    return Model([input1,input2,input3,input4,input5,input6],outputs)

def CDR_model_original3(mlp_units,n_classes,rate=0.1,max_length=30, multi_class=False):
    input1=Input(shape=(30,))
    input2=Input(shape=(30,))
    input3=Input(shape=(30,))

    encoder1=build_t_encoder(num_layers=4,d_model=256,num_heads=4,dff=512,input_vocab_size=22,maximum_position_encoding=max_length)(input1)
    encoder2=build_t_encoder(num_layers=4,d_model=256,num_heads=4,dff=512,input_vocab_size=22,maximum_position_encoding=max_length)(input2)
    encoder3=build_t_encoder(num_layers=4,d_model=256,num_heads=4,dff=512,input_vocab_size=22,maximum_position_encoding=max_length)(input3)

    x = Concatenate()([encoder1, encoder2,encoder3])
    for dim in mlp_units:
        x = Dense(dim, activation="relu")(x)
        x = Dropout(rate)(x)
    act = 'softmax' if multi_class else 'sigmoid'
    outputs = Dense(n_classes, activation=act)(x)
    return Model([input1,input2,input3],outputs)

def CDR_model_single(mlp_units=[512,128,64],n_classes=1,rate=0.1,max_length=30):
    shape = (max_length,)
    input1 = Input(shape=shape)

    encoder1 = build_t_encoder(num_layers=4,d_model=256,num_heads=4,dff=512,input_vocab_size=22,maximum_position_encoding=max_length,shape=shape)(input1)

    x = encoder1

    for dim in mlp_units:
        x = Dense(dim, activation="relu")(x)
        x = Dropout(rate)(x)
    act = 'softmax' if n_classes>1 else 'sigmoid'
    outputs = Dense(n_classes, activation=act)(x)
    
    return Model(input1,outputs)

def CDR_model_dense(mlp_units=[512,128,64],n_classes=1,rate=0.1,max_length=30):
    shape = (max_length,)
    input1=Input(shape=shape)

    x = Dense(256)(input1)
    x = Dense(512)(x)
    x = Dense(512)(x)

    for dim in mlp_units:
        x = Dense(dim, activation="relu")(x)
        x = Dropout(rate)(x)
    act = 'softmax' if n_classes>1 else 'sigmoid'
    outputs = Dense(n_classes, activation=act)(x)

    return Model(input1,outputs)