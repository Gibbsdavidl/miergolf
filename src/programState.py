


def initProgramState (args):
    configText  = open(args[0], 'r').read().split("\n")
    programArgs = filter (lambda zed: zed is not None , map(lambda x: parser(x), configText)) 
    programState = {"runs" : int(programArgs[0]),
                    "ants" : int(programArgs[1]),
                    "ct"   : float(programArgs[2]),
                    "local": int(programArgs[3]),    # perform local optimization?  1 == yes, 0 == no
                    "evap" : float(programArgs[4]),  # evaporation rate on pheromone
                    "damp" : float(programArgs[5]),  # probabilities dampening
                    "alph" : float(programArgs[6]),  # put weight on the edge weights
                    "beta" : float(programArgs[7]),  # put weight on the pheromone weights
                    "mode" : programArgs[8],         # tx (transmit), rx (receive), both
                    "opton": programArgs[9],         # optimize on
                    "tx"   : float(programArgs[10]), # tx cutoff for counting if a node is "touched"
                    "rx"   : float(programArgs[11]), # rx cutoff for counting if a node is "touched"
                    "bestEver" : (0.0, 0.0, []),
                    "bestRest" : (0.0, 0.0, []),
                    "bestIter" : (0.0, 0.0, []),
                    "config" : args[0],              # the config file name
                    "graphfile" : args[1],           # the graph file name 
                    "c" : 1.0,                       # current distance to convergence
                    "k"    : float(programArgs[12]), # size of the set we're finding
                    "cpus"  : int(programArgs[13]),  # number of cpus to use 
                    "lineGraph" : programArgs[14],   # some day will perform this on regular graphs too
                    "permutes" : int(programArgs[15]), # for running the permutation mode
                    "ptoKratio" : 0.0,                 
                    "runsDone" : 0,
                    "iters" : 0}
    return(programState)


def initSimState (args, s):
    configText  = open(args[0], 'r').read().split("\n")
    programArgs = filter (lambda zed: zed is not None , map(lambda x: parser(x), configText)) 
    s["nodes"] = int(programArgs[16])
    s["edges"] = int(programArgs[17])
    s["degpow"] = int(programArgs[18]) # the degree power of the power law graph
    s["timesteps"] = int(programArgs[19])
    s["disp"] = float(programArgs[20])
    s["full"] = programArgs[21]  # should be fullrun, which does ants, or anything else, and just the sim is done
    s["storesize"] = programArgs[22] # the size of the node stores
    return(s)


def parser(x):
    a = x.strip()
    b = a.split(':')
    if len(b) > 1:
        return (b[1].strip())

