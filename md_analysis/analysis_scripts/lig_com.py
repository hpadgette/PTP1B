#!/ usr / bin / env python

# author: Anika J. Friedman

import mdtraj as md
import argparse
import math
import random
import numpy as np
import sys

#Declare arguments
parser = argparse.ArgumentParser(description = 'Concatenate uncorrelated frames from multiple trajectories')
parser.add_argument('-d', required=True, type=str, help='Directory Path for PTP1B repository')
parser.add_argument('-l', required=True, type=str, help='Ligand Name')
parser.add_argument('-c', required=False, default=False, type=bool, help='Are multiple ligands present?')

#Import Arguments
args = parser.parse_args()
directory = args.d
lig = args.l
combo = args.c

#Import custom modules
sys.path.insert(1, directory + '/util/')
import mdfunc

#Load trajectory file names
if combo == True:
    traj_gro = open('traj_gro_' + lig + '_combo_file.txt', 'r').readlines()[0].split(' ')
    traj_xtc = open('traj_xtc_' + lig + '_combo_file.txt', 'r').readlines()[0].split(' ')
    uncorr = open('uncorr_file_' + lig + '_combo.txt', 'r').readlines()[0].split(' ')
else:
    traj_gro = open('traj_gro_' + lig + '_file.txt', 'r').readlines()[0].split(' ')
    traj_xtc = open('traj_xtc_' + lig + '_file.txt', 'r').readlines()[0].split(' ')
    uncorr = open('uncorr_file_' + lig + '.txt', 'r').readlines()[0].split(' ')

#Establish empty vector for all RMSD values
dist_all = []

#Load centroid
if combo == True:
    ref_pdb = md.load('Full_' + lig + '_combo_cluster.pdb')
else:
    ref_pdb = md.load('Full_' + lig + '_cluster.pdb')
ref_top = ref_pdb.topology
ref_ns = ref_pdb.atom_slice(ref_top.select('backbone or resname AD or resname BBR')) #Select only atoms in the protein or ligand
ref_top = ref_ns.topology

#Load and concatenate trajectories
for i in range(len(traj_gro)):
    #Load files
    File_gro = traj_gro[i].strip()
    File_traj = traj_xtc[i].strip()
    File_uncorr = uncorr[i].strip()

    #Load trajectory
    traj = md.load(File_traj, top=File_gro)
    top = traj.topology
    if i != 2:
        traj_ns = traj.atom_slice(top.select('backbone or resname AD or resname BBR')) #Select only atoms in the protein or ligand
    else:
        traj_ns = traj.atom_slice(top.select('(backbone or resname AD or resname BBR) and not resid 0')) #Select only atoms in the protein or ligand
    top = traj_ns.topology

    #Limit trajectory to uncorrelated frames
    uncorr_ind_string = open(File_uncorr, 'r').readlines()
    uncorr_ind = np.zeros(len(uncorr_ind_string), dtype=int)
    for j in range(len(uncorr_ind_string)):
        uncorr_ind[j] = int(j)
    traj_uncorr = traj_ns.slice(uncorr_ind)
#    traj_uncorr = traj_ns

    #Align to reference with a3 and a6 helices
    res_align = top.select('resid 185 to 199 or resid 264 to 280')
    
    traj_ns_align = traj_uncorr.superpose(ref_ns, atom_indices = res_align)

    #seperate ligand carbon atoms
    lig_only_ref = ref_ns.atom_slice(ref_top.select('resname ' + str(lig))) #reference
    lig_only_traj = traj_ns_align.atom_slice(top.select('resname ' + str(lig))) #trajectory

    lig_only_ref_top = lig_only_ref.topology
    lig_only_traj_top = lig_only_traj.topology
        
    #Compute COM of ligand
    com = md.compute_center_of_mass(lig_only_traj)
    com_ref = md.compute_center_of_mass(lig_only_ref)
        
    #Compute displacment
    time, dim = np.shape(com)
    displacment = np.zeros(time)
    for j in range(time):
        displacment[j] = (com[j][0] - com_ref[0][0])**2 + (com[j][1] - com_ref[0][1])**2 + (com[j][2] - com_ref[0][2])**2
    print(math.sqrt(np.mean(displacment)))

    #Add to full RMSD list
    for j in displacment:
        dist_all.append(j)

    print('Trajectory ' + str(i) + ' loaded')

#Random sample the full RMSD 20x
sample = 20
RMSD_bootstrap = np.zeros(sample)
for i in range(sample):
    dis_sample = random.sample(dist_all, 200)
    RMSD_bootstrap[i] = math.sqrt(np.mean(dis_sample))

print('mean: ' + str(np.mean(RMSD_bootstrap)))
print('max: ' + str(max(RMSD_bootstrap)))
print('min: ' + str(min(RMSD_bootstrap)))

if combo == True:
    np.savetxt('RMSD_bootstrap_' + lig + '_combo.txt', RMSD_bootstrap)
else:
    np.savetxt('RMSD_bootstrap_' + lig + '.txt', RMSD_bootstrap)

