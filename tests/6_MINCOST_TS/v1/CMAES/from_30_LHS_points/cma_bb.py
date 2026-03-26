#!/usr/bin/env python
import sys
import os
import numpy as np  # for median, zeros, random, asarray
from numpy import loadtxt
import subprocess


### solver imports (add other imports if necessary)
try: import cma
except: pass  # may not be installed

fmin = cma.fmin_con2
n=5
m=6

bbe=0
last_x=n*[0.0]
obj=1e20
cstr=m*[0.0]


def c(x):
    global id, bbe, last_x, cstr, obj
    last_x=x
    bbe=bbe+1
    obj=1e20

    #print("Constraints: bbe=",bbe," x=",x)
    #print("Constraints: bbe=",bbe," bbx=",bbx)

        
    # blackbox evaluation:
    XfileName = 'x_tmp_'+id+'.txt'
    SolarOutFileName = 'solar_output_tmp_'+id+'.txt'

    with open(XfileName, 'w') as x_file:
        for i in range(n):
            print(x[i]," ",file=x_file,end="")
        x_file.close()

    os.system("./solar 6 "+XfileName+" > "+SolarOutFileName)

    with open(SolarOutFileName, 'r') as solar_output:
        s=solar_output.read()
        solar_output.close()

    s.strip("\n")
    tmp=s.split(" ")

    obj=float(tmp[0])
    
    for j in range(m):
        cstr[j]=float(tmp[j+1])
   
    OutFileName = 'stats_'+id+'.txt'

    if bbe==1:
        with open(OutFileName, 'a') as out:
            print("bbe bbo sol",file=out)
            out.close()

    with open(OutFileName, 'a') as out:
        print(bbe," ",tmp," ",x,file=out)
        out.close()

    # print(tmp)

    return cstr

def f(x):
    global id, bbe, last_x, cstr, obj

    for i in range(n):
        if last_x[i] != x[i]:
            raise Exception("Error: not the right x in bb.f(): ",x," except of ",last_x)
        
    # print('Obj:',obj)

    return obj



id = sys.argv[1]

X0FileName = "X0_"+id+".txt"
line = loadtxt(X0FileName, unpack=False, max_rows=1)
x0 = line.tolist()

if len(x0) != n:
    raise Exception("Error: size if x0 not equal to n. X0:",x0," N:", n)

lb =  [793.0,    2.0,    2.0,    0.01,     0.01]
ub =  [995.0,   50.0,   30.0,    5.00,     5.00] 

fmin(f, x0, 1.5, constraints=c,find_feasible_first=True, find_feasible_final=True,
    options = {"bounds": [lb, ub],
                "maxfevals": str(3000),
                "seed": 1,
                 "verbose": -9},
    eval_initial_x=True,restarts=0)
        
