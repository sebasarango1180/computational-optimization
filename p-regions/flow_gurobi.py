# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 13:49:25 2019

@author: matrio
"""

import numpy as np
import os
from gurobipy import *
import math
import time

from data import order_pot
from data import poten
from data import readdata

# %%


def MathModel(c, d, n, p, mod_set, ii, Ic):

    ct=time.time()
    
#    Parameters definition

    set_I = range(n)
    set_K = range(p)
    set_N = []
    for i in range(n):
        Ni = []
        for j in range(n):
            if c[i][j] == 1:
                Ni.append(j)
        set_N.append(Ni.copy())

    # Modelo flow01_instance_model
    mod = Model(f"flow01_{ii}_{Ic}")


# %%
#    Variables definition.
    vart = mod.addVars(n, n, vtype=GRB.BINARY, name='t')  # t_ij
    vary = mod.addVars(n, p, vtype=GRB.BINARY, name='y')  # y_ik
    varw = mod.addVars(n, p, vtype=GRB.BINARY, name='w')  # w_ik
    varf = mod.addVars(n, n, p, lb = 0, vtype=GRB.CONTINUOUS,
                       name='f')  # f_ijk

    mod.update()

# %%
##### ORIGINAL CONSTRAINTS #####


    # Constraint set 1 eq 42
    for i in set_I:
        eq42 = LinExpr()
        for k in set_K:
            eq42.addTerms(1, vary[i, k])
        mod.addConstr(eq42 == 1)


    # Constraint set 2 eq 43
    for i in set_I:
        for k in set_K:
            mod.addConstr(varw[i,k] <= vary[i, k])
            

    # Constraint set 3 eq 44
    for k in set_K:
        eq_44 = LinExpr()
        for i in set_I:
            eq_44.addTerms(1,varw[i, k])
        mod.addConstr(eq_44==1)


    # Constraint sets 4 and 5, eq 45 and 46
    for i in set_I:
        for j in set_N[i]:
            for k in set_K:
                mod.addConstr(varf[i,j,k] <= (n - p) * vary[i,k])
                mod.addConstr(varf[i,j,k] <= (n - p) * vary[j,k])


    # Constraint set 6 eq 47
    for i in set_I:
        for k in set_K:
            sum47_1 = LinExpr()
            for j in set_N[i]:
                sum47_1.addTerms(1, varf[i, j, k])

            sum47_2 = LinExpr()
            for j in set_N[i]:
                sum47_2.addTerms(1, varf[j, i, k])
            
            mod.addConstr(sum47_1 - sum47_2 >= vary[i, k] - (n - p) * varw[i, k])

    # Constraint set 7 eq 48
    for i in set_I:
        for j in set_I:
            if i <= j:
                for k in set_K:
                    mod.addConstr(vart[i, j] >= vary[i, k] + vary[j, k] - 1)


# %%
##### OBJECTIVE FUNCTION #####

#    Objective function eq 1
    mod.setObjective(quicksum(d[i][j]*vart[i, j]
                              for i in set_I for j in set_I if j > i),
                     GRB.MINIMIZE)

    # Solve
    mod.setParam(GRB.Param.OutputFlag, 0)
    mod.setParam(GRB.Param.TimeLimit, 30)

    mod.update()

    mod.write('Flow_' + str(Ic) + '.mps')
#    print(str(ii) + "\t" + str(Ic))
#
#    mod.optimize()
#
#    ct=time.time()-ct
#    print(mod.status,"\t",mod.objVal,"\t",ct)

#    print()
#    for variable in mod.getVars():
#        if variable.X > 0.1:
#            print(variable.getAttr('varName') + "\t" + str(variable.X))
#
#    print()
#    print(mod.status,"\t",mod.objVal,"\t",ct)

    return mod


# %%
# Number of instances
Is = 10

# Number of models, where model 0 is the non altered basic order model
M = 1

# Combinations
Combs = order_pot(poten(list(range(1, M))))

parN,parP,parC,parD=readdata()

Ic = 0
h=0

# For instance ii in set of instances Is
for ii in range(Is):
    # For model set mn in set of models:
    for mn in Combs:
        Ic = Ic + 1
        h = h + 1
        MathModel(parC[ii], parD[ii], parN[ii], parP[ii], mn, ii, Ic)
        print(h,"\t",ii,"\t",mn)
        print("")
