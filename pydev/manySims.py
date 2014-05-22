'''
number of restarts   : 3
number of ants       : 20
converge threshold   : 0.0001
local optimization   : 21
evaporation rate     : 0.1
dampening            : 0.99
alpha                : 1.0
beta                 : 1.0
mode                 : both
optimize on          : touch
transmit threshold   : 0.000001
receive threshold    : 0.000001
k                    : 4
cpus                 : 4
lineGraph            : linegraph
permutes             : 1000
sim nodes            : 15
sim edges            : 15
sim degree power     : 6
timesteps            : 1000
disipate             : 0.1
full run             : fullrun
'''

def doMainSims():
    import subprocess
    numRestarts="3"
    numAnts="20"
    converge="0.0001"
    local="100"
    evap="0.1"
    damp="0.99"
    alpha="1.0"
    beta="1.0"
    mode="both"
    opton="touch"
    tx="0.000001"
    rx="0.000001"
    k="4"
    cpus="4"
    linegraph="linegraph"
    perms="1000"
    simnode="100"
    simedges="100"
    simdeg="6"
    timesteps="5000"
    dissipate="0.2"
    full="fullrun"

    

    for simi in xrange(25):
        print simi
        graphname = "simgraph" + str(simi) + ".txt"
        outname = "output_" + str(simi) + ".txt"
        theCmd = "ipython ../mipdao/pydev/flowSim.py configSim.txt " + graphname + " > " + outname
        subprocess.call(theCmd, shell=True)
