#!/usr/local/bin/python3
import argparse
import numpy as np
import seaborn as sns
import pandas as pd
import sys
import os
# average over window, window has to be odd
def window_avg(values, window_size):

    # check for odd window
    if window_size % 2 == 0:
        sys.error("window length has to be odd number")

    half = int((window_size-1)/2)

    # values are a list of floats
    lengthened = values.copy()

    # add half a window to beginning with value of the first element
    for i in range(0, half):
        lengthened.insert(0, values[0])

    # add half a window to beginning with value of the last element
    for i in range(0, half):
        lengthened.append(values[-1])

    avg = []
    # sum values over list, creating a new list
    for i in range(0, len(values)):

        this_sum = 0
        for j in range(0, window_size):
            this_sum += lengthened[i+j]
        avg.append(this_sum / window_size)

    return avg

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--saliency', '-i', type=str, default=["315-04-1D02", "56.a.09_Heavy", "MEDI8852-V", "315-02-1F07", "39.29NCv1", "PN-SIA28", "429B01", "FI6v3", "CT149", "S9-3-37", "CR9114", "F10", "CR6261", "27F3", "315-27-1C08", "CR8020", "16.g.07_Heavy", "315-53-1A09", "31.a.83_Heavy", "315-13-1B02", "CR8043", "16.a.26_Heavy", "31.b.09_Heavy"],nargs='+',help='JSON file with saliency maps.')
    parser.add_argument('--pdb', '-pdb', type=str, default=["5WCD", "5K9K", "5JW3", "5WCC", "4KVN", "8GV6", "6NZ7", "3ZTJ", "4R8W", "6E3H", "4FQI", "3FKU", "3GBN", "5WKO", "5WCA", "3SDY", "5KAN", "5U4R", "5KAQ", "5TY6", "4NM8", "5K9Q", "5K9O"],nargs='+', help='PDB id to fetech')
    parser.add_argument('--chain', '-c', type=str, default=['H','H','G','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H'],nargs='+', help='PDB chain')
    parser.add_argument('--window', '-w', type=int, default=3, help='Sliding window length. Default 3.')
    parser.add_argument('--saliency_path','-p', type=str, default='../result/average_saliency/',help="input saliency tsv file dir")
    parser.add_argument('--output_path','-o', type=str, default='graph/pdb_saliency/',help="output saliency png")
    args = parser.parse_args()
    # os.makedirs(args.output_path, exist_ok=True)
    window = args.window
    with open('pymol_viz_test.py', 'w') as f:
        f.write("#!/usr/bin/env python\n")
        f.write("import time\n")
        f.write("from pymol import stored\n")
        f.write("import pymol\n\n")
        f.write("surface_color = 'white'\n")
        f.write("antibody_colors = 'silver_yellow_red'\n")
        f.write("pymol.cmd.bg_color( \'white\' )\n")
        f.write("pymol.cmd.viewport( \'2000\', \'2000\' )\n")
        f.write("pymol.cmd.set( \'sphere_scale\', \'0.5\' )\n")
        for i,pdb in enumerate(args.pdb):
            salieny_f = f'{args.saliency_path}{args.saliency[i]}_gradcam.tsv'
            gradcam_values = pd.read_csv(salieny_f,sep='\t').loc[:,args.saliency[i]].tolist()[:120]
            new_gradcam = window_avg(gradcam_values, window)
        

            # write script to run in pymol
            f.write("pymol.cmd.fetch( \'" + pdb + "\' )\n")
            f.write("time.sleep( 1 )\n")
            f.write("pymol.cmd.hide( 'all')\n")
            f.write("pymol.cmd.show( 'cartoon','chain H + chain L' )\n")
            f.write("pymol.cmd.show( 'surface','chain A + chain B' )\n")
            f.write("pymol.cmd.set( 'transparency', '0.5')\n")
            f.write("pymol.cmd.set( 'surface_color', surface_color)\n")
            f.write("pymol.cmd.set( 'surface_mode' , '3')\n")
            f.write("pymol.cmd.set( 'ray_transparency_oblique')\n")
            f.write("pymol.cmd.set( 'transparency_mode', '1')\n")
            f.write("pymol.cmd.set( 'sphere_scale', '1' )\n")
            f.write(f"seq = \'\'.join(pymol.cmd.get_fastastr(\'chain {args.chain[i]}\').split(\'\\n\')[1:])\n")
            # alter b-factor columns with custom coloring
            cam_string = str(", ".join(map(str, new_gradcam)))
            f.write(f"pymol.cmd.alter(\'{pdb}\', \'b=0.0\' )\n")
            f.write(f"stored.cam = [ {cam_string} ]\n")
            f.write("stored.cam += [0.0] * (len(seq) - len(stored.cam))\n")
            f.write(f"pymol.cmd.alter( \'chain {args.chain[i]} and n. CA\', \'b=stored.cam.pop(0)\' )\n")

                # a bunch of nice color schemes
                # f.write( "pymol.cmd.spectrum(\'b\', \'slate_yellow_red\', \'" + lchains[c] + " and n. CA\' )\n")
                # f.write( "pymol.cmd.spectrum(\'b\', \'silver_yellow_red\', \'" + lchains[c] + " and n. CA\' )\n")
                # f.write( "pymol.cmd.spectrum(\'b\', \'palecyan_silver_magenta\', \'" + lchains[c] + " and n. CA\' )\n")
                # f.write( "pymol.cmd.spectrum(\'b\', \'aquamarine_silver_red\', \'" + lchains[c] + " and n. CA\' )\n")
                # f.write( "pymol.cmd.spectrum(\'b\', \'lightblue_silver_red\', \'" + lchains[c] + " and n. CA\' )\n")
            f.write(f"pymol.cmd.spectrum( \'b\', antibody_colors, \'{pdb} and n. CA\' )\n")

                # ligands, water, etc
            f.write(f"pymol.cmd.select( \'ligs_{pdb}\', \'{pdb} and het\' )\n")
            f.write("pymol.cmd.select( \'water\', \'resn hoh\' )\n")
            f.write("pymol.cmd.remove( \'water\' )\n")
            f.write(f"pymol.cmd.remove( \'ligs_{pdb}\' )\n")
                    
            f.write(f"pymol.cmd.set_view(())\n")
            f.write(f"pymol.cmd.png( \'{args.output_path}{pdb}.png"+"\', dpi=300 )\n")
            f.write("pymol.cmd.delete( \'all\')\n")
            f.write("\n")
        f.close()