module LocOpt (makeSolns, chunk, connectA, connects2Dedge, connectB) where

import System.Random
import ProgramState
import Graph
import Utilities

--localOptimize :: State -> Digraph -> ScoreData -> Ant -> Ant
-- we want alternate solutions
--localOptimize s digraph scoredata (node, soln, score)  
 --             = locop s digraph scoredata (node, soln, score) (makeSolns soln digraph)
              

--locop :: State -> Digraph -> ScoreData -> Ant -> [Digraph] -> Ant
-- sort through all the potential variations --  keep the best
--locop s dig scoredata ant [] = ant
--locop s dig scoredata ant (soln:solns) = locop s dig scoredata ant' solns
 --     where ant' = checkBest s ant (updateAntSoln ant soln) scoredata

--updateAntSoln :: Ant -> Digraph -> Ant
--updateAntSoln (n1, es, scorex) soln = (n1, soln, 0.0)


-- for each list in starts -> x
--    for each Dedge in x
--        find the connector to match


midchunk (a,b,c) = (b!!0)


makeSolns :: Digraph -> Digraph -> [Digraph]
--- make lots of potential solutions given a soln (asoln)
makeSolns asoln digraph = solutionList  
          where chunks = chunk [] asoln  :: [([Dedge], [Dedge], [Dedge])]
                starts = map (\(xs, chunk, ys) -> connectA (chunk !! 1) digraph) chunks :: [Digraph]
                zipstarts = zip chunks starts
                connectors = map (\(chunker,startlist) ->
                                 (map (\starter -> connectB starter digraph (midchunk chunker)) startlist)) zipstarts :: [Digraph]
                solutionList = joinSolutions chunks starts connectors ([] :: [Digraph]) :: [Digraph]
                  

joinSolutions:: [([Dedge], [Dedge], [Dedge])] -> [Digraph] -> [Digraph] -> [Digraph] -> [Digraph] 
joinSolutions _ [] [] newSolns = newSolns
joinSolutions ((beg,ch,end):chunks) (st:starts) (con:cons) newSolns = joinSolutions chunks starts cons grownSolns
              where grownSolns = (processStarts beg end con st ([]::[Digraph])) ++ newSolns

-- NEED ONE MORE LEVEL HERE TO PROCESS THE LIST OF STARTERS FOR EACH CHUNK... --
processStarts :: Digraph -> Digraph -> Digraph -> Digraph -> [Digraph] -> [Digraph]
processStarts _ _ [] [] l = l
processStarts begin end (con:connector) (str:starts) l = 
              processStarts begin end connector starts ((begin ++ (con:str:end)) : l)

chunk :: [a] -> [a] -> [([a], [a], [a])]
chunk [] []           = [([],   [],   [])]
chunk [] (y:x:[])     = [([],   [y,x],[])]
chunk begin (x:y:[])  = [(begin,[x,y],[])]
chunk [] (x:y:end)    = ([], [x,y], end)   : (chunk [x] (y:end))
chunk begin (x:y:end) = (begin, [x,y],end) : (chunk (begin++[x]) (y:end))


connectA :: Dedge -> Digraph -> Digraph
-- take an edge, A -> B, and return a list of new Dedges, with new destinations A -> X
connectA a dig = filter (\x -> and [(connects2Dedge a x), (notEndingAtStart2 x)]) dig 

connects2Dedge :: Dedge -> Dedge -> Bool   -- Dedges: (going-to)-(coming-from)
connects2Dedge (e1,e2,e3,e4) (a,b,c,d)
         | e2 == b = True
         | otherwise = False

notEndingAtStart2 :: Dedge -> Bool -- Dedge (going-to)-(coming-from)
-- don't go back to the start until the end.  
notEndingAtStart2 (a1,b1,c1,d1) 
  | a1 == (("Start1","Start2",1.0,1.0):: Edge) = False
  | otherwise = True

connectB :: Dedge -> Digraph -> Dedge -> Dedge
-- given the new connector edge, x, what connects to it?
--       starter (X-a)       digraph           chunk *(B - _)*  <- (_ - A)
connectB  (b1,b2,b3,b4) ((d1,d2,d3,d4):dig) (x1,x2,x3,x4) 
         |  d1 == x1 && d2 == b1 = (d1,d2,d3,d4)
         | otherwise = connectB (b1,b2,b3,b4) dig (x1,x2,x3,x4)

