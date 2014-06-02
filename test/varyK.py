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



def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)

    pathtomipdao = args[0]
    graphname = args[1]
    outputpath = args[2]
    cpus = int(args[3])

    for k in [1,2,3,4,5,6,7,8]:
        for i in xrange(3):
            print "K: " + str(k) + "   " + str(i)
            outname = outputpath + "_varyK_output_" + str(k) + "_" + str(i) + ".txt"
            configname = outputpath + "_varyK_config_"+ str(k) + "_" + str(i) + ".txt"
            printConfig(configname, k, cpus, 0)
            theCmd = "ipython " + pathtomipdao + "pydev/main.py "+configname+ " " + graphname + " > " + outname
            os.system(theCmd)



def printConfig(filename, k, cpus, nodes):
    fout = open(filename,'w')
    steps = nodes * 2000
    fout.write("number of restarts   : 3\n")
    fout.write("number of ants       : 18\n")
    fout.write("converge threshold   : 0.0001\n")
    fout.write("local optimization   : 20\n")
    fout.write("evaporation rate     : 0.1\n")
    fout.write("dampening            : 0.99\n")
    fout.write("alpha                : 1.0\n")
    fout.write("beta                 : 1.0\n")
    fout.write("mode                 : both\n")
    fout.write("optimize on          : score\n")
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
    fout.write("dissipate             : 0.01\n")
    fout.write("full run             : fullrun\n")
    fout.write("store size         : 20\n")
    fout.close()



if __name__ == "__main__":
    import sys
    import os
    import getopt
    main()
