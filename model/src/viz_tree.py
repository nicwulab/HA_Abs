import pandas as pd
import tensorflow_decision_forests as tfdf
from ast import literal_eval
from tqdm import tqdm
from collections import Counter
import json

import seaborn as sns
import matplotlib.pylab as plt
import logomaker
from matplotlib import pyplot as plt

def get_rule(node):
    
    char_index = node['condition']['attribute']
    char_list = node['condition']['mask']

    try:
        class_distribution = node['value']['distribution']

        return {'char_index': char_index, 
                'char_list': char_list,
              'class_distribution': class_distribution}
    except KeyError:
        print('KeyError', node)

    return {'char_index': char_index, 
                'char_list': char_list}
    
def depth_traverse(node, node_dict, node_id=(0,0)):
    '''
    recursive depth first search traversal
    node_list: [traversed_nodes]
    node_id: (level_index, children_index)
    '''

    if 'condition' in node:
        if node_id in node_dict:
            level_id, children_id = node_id
            node_id = (level_id, children_id+1)
        node_dict[node_id] = get_rule(node)

        if 'children' in node:
            level_id, children_id = node_id
            for i, c in enumerate(node['children']):
                depth_traverse(c, node_dict, (level_id+1, children_id+i))
    
def convert_tree(node_dict):
    node_dict_tree = {}
    for (node_id, node) in node_dict.items():
        level_id, children_id = node_id
        if level_id not in node_dict_tree:
            node_dict_tree[level_id] = {}
        node_dict_tree[level_id][children_id] = node
    return node_dict_tree

def read_tree(html):
    start_string = 'display_tree({}, '
    end_string = ', "#tree_plot'
    start_pos = html.index(start_string) + len(start_string)
    end_pos = html.index(end_string)
    tree_dict = literal_eval(html[start_pos:end_pos])

    node_dict = {}
    depth_traverse(tree_dict, node_dict)
    return convert_tree(node_dict)

def update_dict_counter(dic, key, cntr):

    if key in dic:
        dic[key].update(cntr)
    else:
        dic[key] = cntr
        
    return 0

def count_occur(node, level_weight):
    
    cntr = Counter(node['char_list'])
    
    class_dist = node['class_distribution']
    # entropy = max(class_dist) / min(class_dist)
    entropy = max(class_dist) / (1 - max(class_dist))
    
    for char, cnt in cntr.items():
        cntr[char] = cnt * entropy * level_weight
        
    return cntr

def update_rule_dict(tree, rule_dict=None):
    
    if rule_dict is None:
        rule_dict = {}
    
    for lv_id, lv in tree.items():
        
        total_level = len(tree.keys())
        level_weight = (total_level + 1 - int(lv_id)) / total_level
        
        for node_id, node in lv.items():
            status = update_dict_counter(rule_dict, node['char_index'], count_occur(node, level_weight))

    return rule_dict

def get_rules_dict(json, rule_dict=None):
    if rule_dict is None:
        rule_dict = {}
    
    for tree_id, tree in json.items():
        rule_dict = update_rule_dict(tree, rule_dict)
            
    return rule_dict

def rule_dict_to_df(rule_dict):
    rule_df = {}
    cdr_char = list('XEDRKHQNSTPGCAVILMFYW')
    cdr_char.append('<OOD>')

    for pos, cntr in rule_dict.items():
        cnt_list = []
        for char in cdr_char:
            if char == 'X' or char == '<OOD>':
                continue 
            if char in cntr:
                cnt_list.append(int(cntr[char]))
            else:
                cnt_list.append(0)
                
        rule_df[int(pos)] = cnt_list

    rule_df = pd.DataFrame.from_dict(rule_df, orient='index', columns=list('EDRKHQNSTPGCAVILMFYW'))
    return rule_df.sort_index()

def extract_tree(tree_model, num_trees=300, max_depth=16):
    
    tree_dict = {}

    for idx in tqdm(range(num_trees)):
        html = tfdf.model_plotter.plot_model(tree_model, idx, max_depth)
        tree_dict[idx] = read_tree(html)

    # with open(json_filename, 'w') as f:
    #     json.dump(tree_dict, f)

    # tree_dict = read_json(os.path.join(json_filename))

    rules_dict = get_rules_dict(tree_dict)
    rule_df = rule_dict_to_df(rules_dict)

    return rule_df

def draw_heatmap(rule_df):
    plt.figure(figsize=(60, 6)) 
    ax = sns.heatmap((rule_df/rule_df.max().max()*100).transpose(), yticklabels=True)#, annot=True)
    plt.show()


def logosequence(logoseq, figsize=[200,50]):
    # create Logo object

    fig, ax = plt.subplots(1,1,figsize=figsize)

    crp_logo = logomaker.Logo(logoseq,
                              ax=ax,
                              color_scheme='NajafabadiEtAl2017')

    # style using Logo methods
    crp_logo.style_spines(visible=False)
    crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
    crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

    # # style using Axes methods
    crp_logo.ax.set_ylabel("$-\Delta \Delta G$ (kcal/mol)", labelpad=-1)
    crp_logo.ax.xaxis.set_ticks_position('none')
    crp_logo.ax.xaxis.set_tick_params(pad=-1)
    