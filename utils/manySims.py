'''
A script to automate running many simulations on different sizes of graph.
'''

def main():
    for n in [20,40,60,80,100]:  # number of steps in growing the graph
        for i in xrange(10):     # number of reps
            print str(n) + "   " + str(i)
            graphname = "simgraph_" + str(n) + "_" + str(i) + ".txt"
            outname = "output_" + str(n) + "_" + str(i) + ".txt"
            configname = "config_"+ str(n) + "_" + str(i) + ".txt"
            printConfig(configname, 4, 9, n)
            theCmd = "ipython ../mipdao/pydev/flowSim.py "+configname+ " " + graphname + " > " + outname
            os.system(theCmd)



def printConfig(filename, k, cpus, nodes):
    fout = open(filename,'w')
    steps = nodes * 1000
    fout.write("number of restarts   : 3\n")
    fout.write("number of ants       : 18\n")
    fout.write("converge threshold   : 0.0001\n")
    fout.write("local optimization   : 20\n")
    fout.write("evaporation rate     : 0.1\n")
    fout.write("dampening            : 0.99\n")
    fout.write("alpha                : 1.0\n")
    fout.write("beta                 : 1.0\n")
    fout.write("mode                 : both\n")
    fout.write("optimize on          : touch\n")
    fout.write("transmit threshold   : 0.0001\n")
    fout.write("receive threshold    : 0.0001\n")
    fout.write("k                    : "+ str(k) + "\n")
    fout.write("cpus                 : "+str(cpus)+"\n")
    fout.write("lineGraph            : linegraph\n")
    fout.write("permutes             : 1000\n")
    fout.write("sim nodes            : "+str(nodes)+"\n")
    fout.write("sim edges            : 0\n")
    fout.write("sim degree power     : 0\n")
    fout.write("timesteps            : "+str(steps)+"\n")
    fout.write("dissipate             : 0.1\n")
    fout.write("full run             : fullrun\n")
    fout.write("store size         : 20\n")
    fout.close()



if __name__ == "__main__":
    import sys
    import os
    main()
