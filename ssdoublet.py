import sys
sys.path.append("/Users/sakuma/PycharmProjects/SEQGEN/venv/PyRosetta4.Release.python36.mac.release-213/setup/")
import numpy as np
counter=0
import pyrosetta as prs
prs.init(options="-in:ignore_unrecognized_res true -ignore_zero_occupancy false")

import math
import numpy as np

import angleutils as au
import pose2features_write_aaseq as p2f
def pose2sseabego(pose=None):
    abego_manager = prs.rosetta.core.sequence.ABEGOManager()
    aaseq=prs.rosetta.core.sequence.Sequence(pose).sequence()

    aa_list=[]
    abego_list = []
    #np.array
    phis = []
    psis = []
    omegas = []
    for i in range(1, pose.size() + 1):
        phi = pose.phi(i)
        psi = pose.psi(i)
        omega = pose.omega(i)
        abego_list.append(abego_manager.index2symbol(abego_manager.torsion2index_level1(
            phi, psi, omega)))
        phis.append(phi)
        psis.append(psi)
        omegas.append(omega)
        aa_list.append(aaseq[i-1])
   # print(abego_list)
    phis = np.array(phis)
    psis = np.array(psis)
    omegas = np.array(omegas)

    DSSP = prs.rosetta.protocols.moves.DsspMover()
    DSSP.apply(pose)
    sss = list(pose.secstruct())
    return sss,abego_list,aa_list


#prs.init()

args = sys.argv
infilename=args[1]
outfilename=args[2]
file=open(outfilename,mode="w")
file.write(
    "pdbname pdbnameid abego ss1 ss2 ss1_start ss1_end loop_start loop_end ss2_start ss2_end NLangle CLangle HHangle NNangle CCangle HHdihedral NLdihedral CLdihedral distloop aa_seq abego_seq sse_seq\n")
with open(infilename) as f:

    #unitcount=0
    for line in f:
        pose = prs.Pose()
        pose = prs.pose_from_pdb(line.split("\n")[0])
        pdbname=line.split("\n")[0].split("/")[-1]
        ss_list,abego_list,aa_list=pose2sseabego(pose)
        #print(aa_list)
        p2f.pose2features(pose,ss_list,abego_list,aa_list=aa_list,pdbname=pdbname,fileobject=file)

    file.close()

