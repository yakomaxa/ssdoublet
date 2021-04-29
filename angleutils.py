from math import sqrt
from math import acos
from math import pi
import numpy as np

def bonds2angle(bondA,bondB):

    v=bondA
    u=bondB
    lenv=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
    lenu=sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])
    unit_u = np.array(u) / lenv
    unit_v = np.array(v) / lenu
    value_angle=acos(min(max(np.dot(unit_v,unit_u),-1.0),1.0))*180/pi
    return(value_angle)


def bonds3dihedral(bond12, bond23, bond34):
    ZERO_JUDGE = 0.00000000001
    v21 = bond12
    v32 = bond23
    v43 = bond34

    #c11 = v21[1] * v21[1] + v21[2] * v21[2] + v21[3] * v21[3]
    #c22 = v32[1] * v32[1] + v32[2] * v32[2] + v32[3] * v32[3]
    #c33 = v43[1] * v43[1] + v43[2] * v43[2] + v43[3] * v43[3]

    c11 = v21[0] * v21[0] + v21[1] * v21[1] + v21[2] * v21[2]
    c22 = v32[0] * v32[0] + v32[1] * v32[1] + v32[2] * v32[2]
    c33 = v43[0] * v43[0] + v43[1] * v43[1] + v43[2] * v43[2]


    #c12 = v21[1] * v32[1] + v21[2] * v32[2] + v21[3] * v32[3]
    #c13 = v21[1] * v43[1] + v21[2] * v43[2] + v21[3] * v43[3]
    #c23 = v32[1] * v43[1] + v32[2] * v43[2] + v32[3] * v43[3]

    c12 = v21[0] * v32[0] + v21[1] * v32[1] + v21[2] * v32[2]
    c13 = v21[0] * v43[0] + v21[1] * v43[1] + v21[2] * v43[2]
    c23 = v32[0] * v43[0] + v32[1] * v43[1] + v32[2] * v43[2]

    t1 = c12 * c23 - c13 * c22
    t3 = c11 * c22 - c12 * c12
    t4 = c22 * c33 - c23 * c23

    if (t3 <= ZERO_JUDGE):
        t3 = ZERO_JUDGE
    if (t4 <= ZERO_JUDGE):
        t4 = ZERO_JUDGE

    t3t4 = sqrt(t3 * t4)
    co_dih = t1 / t3t4

    if (co_dih > 1.0):
    # then ! when co > 1 ,or < 1
       co_dih = 1.0
    elif (co_dih < -1.0):
        co_dih = -1.0

    #zahyokei = v21[1] * v32[2] * v43[3] + v32[1] * v43[2] * v21[3] + v43[1] * v21[2] * v32[3] - v43[1] * v32[2] * v21[
     #   3] - v21[1] * v43[2] * v32[3] - v32[1] * v21[2] * v43[3]
    zahyokei = v21[0] * v32[1] * v43[2] + v32[0] * v43[1] * v21[2] + v43[0] * v21[1] * v32[2] - v43[0] * v32[1] * v21[
        2] - v21[0] * v43[1] * v32[2] - v32[0] * v21[1] * v43[2]

    si_dih = sqrt(1.0 - co_dih ** 2)

    if (si_dih < ZERO_JUDGE):
        si_dih = ZERO_JUDGE

    if (zahyokei < 0.0):
        si_dih = - si_dih

    if (zahyokei > 0.0):
        dih_angle = acos(co_dih) * 180 / pi
    else:
        dih_angle = - acos(co_dih) * 180 / pi

    return (dih_angle)

