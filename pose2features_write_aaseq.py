import numpy as np
import angleutils as au
import math
def pose2features(pose=None,ss_list=None,abego_list=None,aa_list=None,pdbname=None,fileobject=None):
    ssflag = 0
    sseStart = 0
    sseEnd = 0
    sseprev = "X"
    abego = []
    sse = []
    dummy = 0

    aaseeq=pose.sequence()


    loopStart = dummy
    loopEnd = dummy
    sseStart = dummy
    sseEnd = dummy
    sseStartEnd = [[sseStart, sseEnd]]
    loopStartEnd = [[loopStart, loopEnd]]
    sstypes = [[sseprev]]
    ssecount = 0
    loopcount = 0
    prev = "L"
    loopabego = [["X", "X", "X"]]
    loopaaseq = [["X", "X", "X"]]

    for resi in range(1, pose.size()):
        ssnow = ss_list[resi]
        ssold = ss_list[resi - 1]
        if (ssold == "L" and ssnow == "L"):
            # loop to loop
            abego.append(abego_list[resi])
        elif (ssold == "L" and ssnow != "L"):
            ssecount = ssecount + 1
            # loop to sse
            sseStart = resi
            loopEnd = resi - 1
            loopStartEnd.append([loopStart, loopEnd])
            #print(abego)
            abego = []
            loopabego.append(abego)
            # sse.append(ss_list[resi])
        elif (ssold != "L" and ssnow == "L"):
            loopcount = loopcount + 1
            # sse to loop
            loopStart = resi
            sseEnd = resi - 1
            sseStartEnd.append([sseStart, sseEnd])
            if (ss_list[sseStart] == ss_list[sseEnd - 1]):
                sstypes.append([ss_list[sseStart]])
            else:
                #print("sse differs")
                sstypes.append(["X"])

            # sseStartEnd.append([loopEnd, sseStart])
            abego.append(abego_list[resi])
            # print(sse[ssecount-1],sse[ssecount])
        # print(resi)

        elif (ssnow != "L" and ssold != "L"):
            # sse to loop
            a = 1
            # print("")
    #print(sseStartEnd)
    count=0
    for ii in range(1, len(loopabego) ):
        count=count+1
        if ("X" in loopabego[ii - 1] or "X" in sstypes[ii - 1]
            or "X" in loopabego[ii] or "X" in sstypes[ii]):
            #print("X contained, skip")
            writeflag=0
            ddummyfa=0
        else:
            if (sstypes[ii - 1] == ["E"] and sstypes[ii] == ["E"]):
                writeflag = 0
                if (ss_list[loopStartEnd[ii][0] - 3] == "E" and
                        ss_list[loopStartEnd[ii][0] - 2] == "E" and
                        ss_list[loopStartEnd[ii][0] - 1] == "E" and
                        ss_list[loopStartEnd[ii][1] + 1] == "E" and
                        ss_list[loopStartEnd[ii][1] + 2] == "E" and
                        ss_list[loopStartEnd[ii][1] + 3] == "E"
                ):
                    writeflag = 1
                    aa_seq = ''.join((aa_list[(loopStartEnd[ii][0] - 3):loopStartEnd[ii][1] + 4]))
                    abego_seq = ''.join(abego_list[(loopStartEnd[ii][0] - 3):loopStartEnd[ii][1] + 4])
                    sse_seq = ''.join(ss_list[(loopStartEnd[ii][0] - 3):loopStartEnd[ii][1] + 4])
                    #print(aa_seq)
                    #print(abego_seq)
                    #print(sse_seq)
                    #print(aa_list[(loopStartEnd[ii][0] - 3):loopStartEnd[ii][1] + 3])
                    # pose.residue().xyz("CA")[0
                    x0_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[0]
                    y0_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[1]
                    z0_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[2]

                    x1_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[0]
                    y1_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[1]
                    z1_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[2]

                    x2_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[0]
                    y2_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[1]
                    z2_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[2]

                    x0_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[0]
                    y0_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[1]
                    z0_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[2]

                    x1_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[0]
                    y1_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[1]
                    z1_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[2]

                    x2_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[0]
                    y2_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[1]
                    z2_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[2]

                    r_n0 = [(x0_n + x1_n) / 2, (y0_n + y1_n) / 2, (z0_n + z1_n) / 2]
                    r_n1 = [(x1_n + x2_n) / 2, (y1_n + y2_n) / 2, (z1_n + z2_n) / 2]

                    r_c0 = [(x0_c + x1_c) / 2, (y0_c + y1_c) / 2, (z0_c + z1_c) / 2]
                    r_c1 = [(x1_c + x2_c) / 2, (y1_c + y2_c) / 2, (z1_c + z2_c) / 2]

                    v_n = np.array(r_n0) - np.array(r_n1)
                    v_c = np.array(r_c0) - np.array(r_c1)

                    vert_n = np.array([x2_n, y2_n, z2_n]) - np.array(r_n1)
                    vert_c = np.array([x2_c, y2_c, z2_c]) - np.array(r_c1)
                    # b=np.array(x1_c) - np.array(r_c1)
                    v_l = np.array([x2_c, y2_c, z2_c]) - np.array([x2_n, y2_n, z2_n])

                    ## define angles and distance
                    NL_angle = au.bonds2angle(v_n, v_l)
                    CL_angle = au.bonds2angle(v_c, -v_l)
                    NC_angle = au.bonds2angle(v_n, v_c)

                    NN_angle = au.bonds2angle(v_n, vert_n)
                    CC_angle = au.bonds2angle(v_c, vert_c)

                    HH_dih = au.bonds3dihedral(v_n, v_l, v_c)
                    NL_dih = au.bonds3dihedral(v_n, v_l, vert_n)
                    CL_dih = au.bonds3dihedral(v_c, -v_l, vert_c)

                    R_NL = math.sqrt(np.dot(v_l, v_l))


                # au.bonds2angle()
            elif (sstypes[ii - 1] == ["E"] and sstypes[ii] == ["H"]):
                writeflag=0
                if (    ss_list[loopStartEnd[ii][0] - 3] == "E" and
                        ss_list[loopStartEnd[ii][0] - 2] == "E" and
                        ss_list[loopStartEnd[ii][0] - 1] == "E" and
                        ss_list[loopStartEnd[ii][1] + 1] == "H" and
                        ss_list[loopStartEnd[ii][1] + 2] == "H" and
                        ss_list[loopStartEnd[ii][1] + 3] == "H" and
                        ss_list[loopStartEnd[ii][1] + 4] == "H" and
                        ss_list[loopStartEnd[ii][1] + 5] == "H"
                ):
                    writeflag = 1
                    aa_seq = ''.join((aa_list[(loopStartEnd[ii][0] - 3):loopStartEnd[ii][1] + 6]))
                    abego_seq = ''.join(abego_list[(loopStartEnd[ii][0] - 3):loopStartEnd[ii][1] + 6])
                    sse_seq = ''.join(ss_list[(loopStartEnd[ii][0] - 3):loopStartEnd[ii][1] + 6])
                    #print(aa_seq)
                    #print(abego_seq)
                    #print(sse_seq)
                    #print(aa_list[(loopStartEnd[ii][0] - 3):loopStartEnd[ii][1] + 5])
                    # Strand
                    #
                    x0_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[0]
                    y0_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[1]
                    z0_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[2]

                    x1_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[0]
                    y1_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[1]
                    z1_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[2]

                    x2_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[0]
                    y2_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[1]
                    z2_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[2]

                    # HElix
                    x0_c = pose.residue(loopStartEnd[ii][1] + 5).xyz("CA")[0]
                    y0_c = pose.residue(loopStartEnd[ii][1] + 5).xyz("CA")[1]
                    z0_c = pose.residue(loopStartEnd[ii][1] + 5).xyz("CA")[2]

                    x1_c = pose.residue(loopStartEnd[ii][1] + 4).xyz("CA")[0]
                    y1_c = pose.residue(loopStartEnd[ii][1] + 4).xyz("CA")[1]
                    z1_c = pose.residue(loopStartEnd[ii][1] + 4).xyz("CA")[2]

                    x2_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[0]
                    y2_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[1]
                    z2_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[2]

                    x3_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[0]
                    y3_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[1]
                    z3_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[2]

                    x4_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[0]
                    y4_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[1]
                    z4_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[2]

                    r_n0 = [(x0_n + x1_n) / 2, (y0_n + y1_n) / 2, (z0_n + z1_n) / 2]
                    r_n1 = [(x1_n + x2_n) / 2, (y1_n + y2_n) / 2, (z1_n + z2_n) / 2]
                    vert_n = np.array([x2_n, y2_n, z2_n]) - np.array(r_n1)

                    r_c0 = np.array([0.74 * x0_c + x1_c + x2_c + 0.74 * x3_c,
                                     0.74 * y0_c + y1_c + y2_c + 0.74 * y3_c,
                                     0.74 * z0_c + z1_c + z2_c + 0.74 * z3_c]) / 3.48
                    r_c1 = np.array([0.74 * x1_c + x2_c + x3_c + 0.74 * x4_c,
                                     0.74 * y1_c + y2_c + y3_c + 0.74 * y4_c,
                                     0.74 * z1_c + z2_c + z3_c + 0.74 * z4_c]) / 3.48
                    vert_c = np.array([x4_c, y4_c, z4_c]) - np.array(r_c1)

                    v_n = np.array(r_n0) - np.array(r_n1)
                    v_c = np.array(r_c0) - np.array(r_c1)
                    v_l = np.array([x4_c, y4_c, z4_c]) - np.array([x2_n, y2_n, z2_n])

                    ## define angles and distance
                    NL_angle = au.bonds2angle(v_n, v_l)
                    CL_angle = au.bonds2angle(v_c, -v_l)
                    NC_angle = au.bonds2angle(v_n, v_c)

                    NN_angle = au.bonds2angle(v_n, vert_n)
                    CC_angle = au.bonds2angle(v_c, vert_c)

                    HH_dih = au.bonds3dihedral(v_n, v_l, v_c)
                    NL_dih = au.bonds3dihedral(v_n, v_l, vert_n)
                    CL_dih = au.bonds3dihedral(v_c, -v_l, vert_c)

                    R_NL = math.sqrt(np.dot(v_l, v_l))

                    #return NL_angle, CL_angle, NC_angle, NN_angle, CC_angle, HH_dih, NL_dih, CL_dih, R_NL

            elif (sstypes[ii - 1] == ["H"] and sstypes[ii] == ["E"]):
                writeflag = 0
                if (    ss_list[loopStartEnd[ii][0] - 5] == "H" and
                        ss_list[loopStartEnd[ii][0] - 4] == "H" and
                        ss_list[loopStartEnd[ii][0] - 3] == "H" and
                        ss_list[loopStartEnd[ii][0] - 2] == "H" and
                        ss_list[loopStartEnd[ii][0] - 1] == "H" and
                        ss_list[loopStartEnd[ii][1] + 1] == "E" and
                        ss_list[loopStartEnd[ii][1] + 2] == "E" and
                        ss_list[loopStartEnd[ii][1] + 3] == "E"
                ):
                    writeflag = 1
                    aa_seq = ''.join((aa_list[(loopStartEnd[ii][0] - 5):loopStartEnd[ii][1] + 4]))
                    abego_seq = ''.join(abego_list[(loopStartEnd[ii][0] - 5):loopStartEnd[ii][1] + 4])
                    sse_seq = ''.join(ss_list[(loopStartEnd[ii][0] - 5):loopStartEnd[ii][1] + 4])
                    #print(aa_seq)
                    #print(abego_seq)
                    #print(sse_seq)

                    ### HELIX
                    x0_n = pose.residue(loopStartEnd[ii][0] - 5).xyz("CA")[0]
                    y0_n = pose.residue(loopStartEnd[ii][0] - 5).xyz("CA")[1]
                    z0_n = pose.residue(loopStartEnd[ii][0] - 5).xyz("CA")[2]

                    x1_n = pose.residue(loopStartEnd[ii][0] - 4).xyz("CA")[0]
                    y1_n = pose.residue(loopStartEnd[ii][0] - 4).xyz("CA")[1]
                    z1_n = pose.residue(loopStartEnd[ii][0] - 4).xyz("CA")[2]

                    x2_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[0]
                    y2_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[1]
                    z2_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[2]

                    x3_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[0]
                    y3_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[1]
                    z3_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[2]

                    x4_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[0]
                    y4_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[1]
                    z4_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[2]

                    # STRAND
                    x0_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[0]
                    y0_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[1]
                    z0_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[2]

                    x1_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[0]
                    y1_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[1]
                    z1_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[2]

                    x2_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[0]
                    y2_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[1]
                    z2_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[2]

                    r_n0 = np.array([0.74 * x0_n + x1_n + x2_n + 0.74 * x3_n,
                                     0.74 * y0_n + y1_n + y2_n + 0.74 * y3_n,
                                     0.74 * z0_n + z1_n + z2_n + 0.74 * z3_n]) / 3.48

                    r_n1 = np.array([0.74 * x1_n + x2_n + x3_n + 0.74 * x4_n,
                                     0.74 * y1_n + y2_n + y3_n + 0.74 * y4_n,
                                     0.74 * z1_n + z2_n + z3_n + 0.74 * z4_n]) / 3.48

                    r_c0 = [(x0_c + x1_c) / 2, (y0_c + y1_c) / 2, (z0_c + z1_c) / 2]
                    r_c1 = [(x1_c + x2_c) / 2, (y1_c + y2_c) / 2, (z1_c + z2_c) / 2]

                    v_n = np.array(r_n0) - np.array(r_n1)
                    v_c = np.array(r_c0) - np.array(r_c1)

                    v_l = np.array([x2_c, y2_c, z2_c]) - np.array([x4_n, y4_n, z4_n])

                    vert_n = np.array([x4_n, y4_n, z4_n]) - np.array(r_n1)
                    vert_c = np.array([x2_c, y2_c, z2_c]) - np.array(r_c1)

                    ## define angles and distance
                    NL_angle = au.bonds2angle(v_n, v_l)
                    CL_angle = au.bonds2angle(v_c, -v_l)
                    NC_angle = au.bonds2angle(v_n, v_c)

                    NN_angle = au.bonds2angle(v_n, vert_n)
                    CC_angle = au.bonds2angle(v_c, vert_c)

                    HH_dih = au.bonds3dihedral(v_n, v_l, v_c)
                    NL_dih = au.bonds3dihedral(v_n, v_l, vert_n)
                    CL_dih = au.bonds3dihedral(v_c, -v_l, vert_c)

                    R_NL = math.sqrt(np.dot(v_l, v_l))
                #print(NL_angle, CL_angle, NC_angle, NN_angle, CC_angle, HH_dih, NL_dih, CL_dih, R_NL)
                #return NL_angle, CL_angle, NC_angle, NN_angle, CC_angle, HH_dih, NL_dih, CL_dih, R_NL
                #############

            elif (sstypes[ii - 1] == ["H"] and sstypes[ii] == ["H"]):
                writeflag = 0
                if (    ss_list[loopStartEnd[ii][0] - 5] == "H" and
                        ss_list[loopStartEnd[ii][0] - 4] == "H" and
                        ss_list[loopStartEnd[ii][0] - 3] == "H" and
                        ss_list[loopStartEnd[ii][0] - 2] == "H" and
                        ss_list[loopStartEnd[ii][0] - 1] == "H" and
                        ss_list[loopStartEnd[ii][1] + 1] == "H" and
                        ss_list[loopStartEnd[ii][1] + 2] == "H" and
                        ss_list[loopStartEnd[ii][1] + 3] == "H" and
                        ss_list[loopStartEnd[ii][1] + 4] == "H" and
                        ss_list[loopStartEnd[ii][1] + 5] == "H"

                ):
                    writeflag = 1
                    aa_seq=''.join((aa_list[(loopStartEnd[ii][0] - 5):loopStartEnd[ii][1] + 6]))
                    abego_seq=''.join(abego_list[(loopStartEnd[ii][0] - 5):loopStartEnd[ii][1] + 6])
                    sse_seq=''.join(ss_list[(loopStartEnd[ii][0] - 5):loopStartEnd[ii][1] + 6])
                    #print(aa_seq)
                    #print(abego_seq)
                    #print(sse_seq)

                # HELIX
                    x0_n = pose.residue(loopStartEnd[ii][0] - 5).xyz("CA")[0]
                    y0_n = pose.residue(loopStartEnd[ii][0] - 5).xyz("CA")[1]
                    z0_n = pose.residue(loopStartEnd[ii][0] - 5).xyz("CA")[2]

                    x1_n = pose.residue(loopStartEnd[ii][0] - 4).xyz("CA")[0]
                    y1_n = pose.residue(loopStartEnd[ii][0] - 4).xyz("CA")[1]
                    z1_n = pose.residue(loopStartEnd[ii][0] - 4).xyz("CA")[2]

                    x2_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[0]
                    y2_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[1]
                    z2_n = pose.residue(loopStartEnd[ii][0] - 3).xyz("CA")[2]

                    x3_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[0]
                    y3_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[1]
                    z3_n = pose.residue(loopStartEnd[ii][0] - 2).xyz("CA")[2]

                    x4_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[0]
                    y4_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[1]
                    z4_n = pose.residue(loopStartEnd[ii][0] - 1).xyz("CA")[2]

                    # HELIX
                    x0_c = pose.residue(loopStartEnd[ii][1] + 5).xyz("CA")[0]
                    y0_c = pose.residue(loopStartEnd[ii][1] + 5).xyz("CA")[1]
                    z0_c = pose.residue(loopStartEnd[ii][1] + 5).xyz("CA")[2]

                    x1_c = pose.residue(loopStartEnd[ii][1] + 4).xyz("CA")[0]
                    y1_c = pose.residue(loopStartEnd[ii][1] + 4).xyz("CA")[1]
                    z1_c = pose.residue(loopStartEnd[ii][1] + 4).xyz("CA")[2]

                    x2_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[0]
                    y2_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[1]
                    z2_c = pose.residue(loopStartEnd[ii][1] + 3).xyz("CA")[2]

                    x3_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[0]
                    y3_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[1]
                    z3_c = pose.residue(loopStartEnd[ii][1] + 2).xyz("CA")[2]

                    x4_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[0]
                    y4_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[1]
                    z4_c = pose.residue(loopStartEnd[ii][1] + 1).xyz("CA")[2]

                    r_n0 = np.array([0.74 * x0_n + x1_n + x2_n + 0.74 * x3_n,
                                     0.74 * y0_n + y1_n + y2_n + 0.74 * y3_n,
                                     0.74 * z0_n + z1_n + z2_n + 0.74 * z3_n]) / 3.48

                    r_n1 = np.array([0.74 * x1_n + x2_n + x3_n + 0.74 * x4_n,
                                     0.74 * y1_n + y2_n + y3_n + 0.74 * y4_n,
                                     0.74 * z1_n + z2_n + z3_n + 0.74 * z4_n]) / 3.48

                    r_c0 = np.array([0.74 * x0_c + x1_c + x2_c + 0.74 * x3_c,
                                     0.74 * y0_c + y1_c + y2_c + 0.74 * y3_c,
                                     0.74 * z0_c + z1_c + z2_c + 0.74 * z3_c]) / 3.48

                    r_c1 = np.array([0.74 * x1_c + x2_c + x3_c + 0.74 * x4_c,
                                     0.74 * y1_c + y2_c + y3_c + 0.74 * y4_c,
                                     0.74 * z1_c + z2_c + z3_c + 0.74 * z4_c]) / 3.48

                    v_n = np.array(r_n0) - np.array(r_n1)
                    v_c = np.array(r_c0) - np.array(r_c1)

                    v_l = np.array([x4_c, y4_c, z4_c]) - np.array([x4_n, y4_n, z4_n])

                    vert_n = np.array([x4_n, y4_n, z4_n]) - np.array(r_n1)
                    vert_c = np.array([x4_c, y4_c, z4_c]) - np.array(r_c1)

                    ## define angles and distance
                    NL_angle = au.bonds2angle(v_n, v_l)
                    CL_angle = au.bonds2angle(v_c, -v_l)
                    NC_angle = au.bonds2angle(v_n, v_c)

                    NN_angle = au.bonds2angle(v_n, vert_n)
                    CC_angle = au.bonds2angle(v_c, vert_c)

                    HH_dih = au.bonds3dihedral(v_n, v_l, v_c)
                    NL_dih = au.bonds3dihedral(v_n, v_l, vert_n)
                    CL_dih = au.bonds3dihedral(v_c, -v_l, vert_c)

                    R_NL = math.sqrt(np.dot(v_l, v_l))
                else:
                    writeflag=0

            if(writeflag==1):
                fileobject.write(pdbname + " " + pdbname + "_" + str(count).zfill(10) + " " +
                                 ''.join(loopabego[ii - 1]) + " " +
                                 ''.join(sstypes[ii - 1]) + " " +
                                 ''.join(sstypes[ii]) + " " +
                                 str(pose.pdb_info().number(sseStartEnd[ii - 1][0])) + " " +
                                 str(pose.pdb_info().number(sseStartEnd[ii - 1][1])) + " " +
                                 str(pose.pdb_info().number(loopStartEnd[ii][0])) + " " +
                                 str(pose.pdb_info().number(loopStartEnd[ii][1])) + " " +
                                 str(pose.pdb_info().number(sseStartEnd[ii][0]))+ " " +
                                 str(pose.pdb_info().number(sseStartEnd[ii][1])) + " " +
                                 str(NL_angle) + " " + str(CL_angle) + " " +
                                 str(NC_angle) + " " + str(NN_angle) + " " + str(CC_angle) + " " +
                                 str(HH_dih) + " " + str(NL_dih) + " " + str(CL_dih) + " " + str(R_NL)+
                                 " " + str(aa_seq) + " "+str(abego_seq)+ " "+str(sse_seq)+"\n")
