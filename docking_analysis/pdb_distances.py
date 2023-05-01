#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 19:35:24 2021

@author: hannahpadgette
"""

# Hannah Padgette
# Shirts Group
# PTP1B Project

# pdb_distances
# input: protein-ligand complex of PTP1B and ligand (designed to be run with alpha-bisabolene)
# output: distance between ligand (residue) and various important residues in the protein (alpha-helices 3-6)

# important required Python packages
import mdtraj as md
import numpy as np
from itertools import product
import argparse

# declare arguments
parser = argparse.ArgumentParser(description = 'Determination of protein-ligand distances between alpha-bisabolene and helices of PTP1B')
parser.add_argument('-p', required=True, help='File name for input structure file (pdb format)')
parser.add_argument('-f', required=True, help='Base for all file names')

# import arguments
args = parser.parse_args()
file_pdb = args.p + '.pdb'
file_base = args.f

# load pdb file
pdb = md.load(file_pdb)

# set residue groups for helices

if pdb.n_residues == 300 or pdb.n_residues == 288: # if pdb starts at res 1
    group_l = [299] # ligand
    group_3 = [186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200] # a3 helix
    group_4 = [221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238] 
    group_5 = [245, 246, 247, 248, 249, 250, 251, 252]
    group_6 = [264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280] # a6 helix
    group_bend = [281, 282, 283, 284, 285, 286]
    group_7 = [287, 288, 289, 290, 291, 292, 293, 294, 295] # a7 helix

if pdb.n_residues == 299 or pdb.n_residues == 287: # if pdb starts at res 0
    group_l = [298] # ligand
    group_3 = [185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199] # a3 helix
    group_4 = [220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237]
    group_5 = [244, 245, 246, 247, 248, 249, 250, 251]
    group_6 = [263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279] # a6 helix
    group_bend = [281, 282, 283, 284, 285, 286]
    group_7 = [286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297] # a7 helix
    
# create residue pairs for distance calculations

pair_a3 = list(product(group_l, group_3))
pair_a4 = list(product(group_l, group_4))
pair_a5 = list(product(group_l, group_5))
pair_a6 = list(product(group_l, group_6))
pair_bend = list(product(group_l, group_bend))

if pdb.n_residues == 300 or pdb.n_residues == 299:
    pair_a7 = list(product(group_l, group_7))

# checking to make sure that residue pairs/groups were created properly 
groups = open('groups.txt', 'w')
groups.write('pair_a3 = ' + str(pair_a3) + '\n\npair_a4 = ' + str(pair_a4) + '\n\npair_a5 = ' + str(pair_a5) + '\n\npair_a6 = ' + str(pair_a6) + '\n\npair_bend = ' + str(pair_bend) + '\n\npair_a7 = ' + str(pair_a7))

# calculate distances


[dist_a3, pairs_a3] = md.compute_contacts(pdb, contacts=pair_a3, scheme='closest', ignore_nonprotein=False, periodic=True, soft_min=False)
[dist_a4, pairs_a4] = md.compute_contacts(pdb, contacts=pair_a4, scheme='closest', ignore_nonprotein=False, periodic=True, soft_min=False)
[dist_a5, pairs_a5] = md.compute_contacts(pdb, contacts=pair_a5, scheme='closest', ignore_nonprotein=False, periodic=True, soft_min=False)
[dist_a6, pairs_a6] = md.compute_contacts(pdb, contacts=pair_a6, scheme='closest', ignore_nonprotein=False, periodic=True, soft_min=False)
[dist_bend, pairs_bend] = md.compute_contacts(pdb, contacts=pair_bend, scheme='closest', ignore_nonprotein=False, periodic=True, soft_min=False)
if pdb.n_residues == 300 or pdb.n_residues == 299:
    [dist_a7, pairs_a7] = md.compute_contacts(pdb, contacts=pair_a7, scheme='closest', ignore_nonprotein=False, periodic=True, soft_min=False)

# calculate average distances
    
avg_dist_a3 = np.average(dist_a3)
avg_dist_a4 = np.average(dist_a4)
avg_dist_a5 = np.average(dist_a5)
avg_dist_a6 = np.average(dist_a6)
avg_dist_bend = np.average(dist_bend)
avg_dist_a7 = np.average(dist_a7)

# write distances to files

dist = open('dist.txt', 'w')
dist.write('dist_a3 = ' + str(dist_a3) + '\n\ndist_a4 = ' + str(dist_a4) + '\n\ndist_a5 = ' + str(dist_a5) + '\n\ndist_a6 = ' + str(dist_a6) + '\n\ndist_bend = ' + str(dist_bend) + '\n\ndist_a7 = ' + str(dist_a7))

dist.write('\n\navg_dist_a3 = ' + str(avg_dist_a3) + '\n\navg_dist_a4 = ' + str(avg_dist_a4) + '\n\navg_dist_a5 = ' + str(avg_dist_a5) + '\n\navg_dist_a6 = ' + str(avg_dist_a6) + '\n\navg_dist_bend = ' + str(avg_dist_bend) + '\n\navg_dist_a7 = ' + str(avg_dist_a7))

print('checkpoint')

