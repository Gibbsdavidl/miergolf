-------------------------------------------------------------------------------------------------------------------------
-- Influence Maximization Problem using hypercube-minmax ant optimization ----------------------------------------------
-- for problems on biological pathway graphs 
-------------------------------------------------------------------------------------------------------------------------

-- main = getArgs >>=  printResult . antopt . buildGraphs . buildState
--                                    ***

-- Input:
-- Tuple with the program state, and the graphs

-- Output:
-- The results.
-------------------------------------------------------------------------------------------------------------------------

module AntOpt (
 optimize
,possibleEdges
,possibleLastEdges
,edgeProb
,cumulate
,getIndex
,maybeSelect
,generateSoln
,localOptimize
,updateDigraph
,convergence
,updateSolns
,edgeupd
,proportional
,inSolutions
,bounded
,edgeElem
,notEndingAtStart
,trimAnt
) where

import ProgramState
import Utilities
import InfoFlow
import Graph
import LocOpt
import System.Random
import Debug.Trace
import Numeric.LinearAlgebra
import Data.Time.Clock.POSIX
import System.IO.Unsafe
import Data.List (nub)

optimize :: (State, Digraph, ScoreData) -> String
optimize (state,digraph,scoredata) 
  | isEmptyState state = "ERROR"  -- stop now!
  | otherwise = showState (minmax state digraph scoredata) scoredata   -- OK, let's compute
 
------------------------------------------------------------------------------------------
-- The particular optimization routine -- MinMax Ant Optimization
-----------------------------------------------------------------------------------------
minmax :: State -> Digraph -> ScoreData -> State
minmax s digraph scoredata
                | doneState s = s  -- all done? , return program state
                | converged s = minmax (newRun s) (resetDigraph digraph) scoredata 
                | otherwise   = minmax s' digraph' scoredata -- or do another iteration --
                                where (digraph', s') = runColony s digraph scoredata

doneState :: State -> Bool
-- are we done with the optimization yet?
doneState s = ((c s) < (ct s)) && ((runsDone s) == (runs s))

converged :: State -> Bool
-- have we converged on a solution?
converged s = ((c s) < (ct s))

newRun :: State -> State
-- going to start a new run... so decrement the runs, and reset the bestreset solution sets
newRun s = s {runsDone = runCount, c = 1, bestRest = (0.0, []), bestIter = (0.0, [])}
           where runCount = ((runsDone s)+1)
      
runColony :: State -> Digraph -> ScoreData -> (Digraph, State)   -- changed states: StdGen,SolnSet,Runs
-- An Ant iteration: run colony needs to return a digraph and a state.
runColony s digraph scoredata = (digraph', s3)
  where (ant, s')  = generateSoln s digraph scoredata       -- generate solutions -- keep the best
        ant2       = localOptimize s' digraph scoredata ant
        s2         = updateSolns ant2 s'      -- is it better than the previous soln?
        digraph'   = updateDigraph digraph s2  -- update the pheromone levels
        s3         = convergence s2 digraph'  -- check for convergence


------------------------------------------------------------------------------------------------
-- compute solutions -- 
-------------------------------------------------------------------------------------------------------------------------

generateSoln :: State -> Digraph -> ScoreData -> (Ant, State)
-- For i = 1 to m ... create an ant ... construct a solution.
--     if the solution is better then a held solution ... keep it.
generateSoln s digraph scoredat = (genSoln (ants s) s digraph scoredat newAnt)

genSoln :: Int -> State -> Digraph -> ScoreData -> Ant -> (Ant,State)
genSoln 0 s digraph scoredata bestAnt = (bestAnt, s)
genSoln currm s digraph scoredata bestAnt = genSoln (currm-1) s' digraph scoredata ant'
          where starnode = ("Start1","Start2",1.0,1.0)::Edge -- start at the START NODE
                newant   = makeNewAnt starnode                    -- our new ant, starting at start.
                (ant,s') = constructSoln ((k s)+1) digraph (newant,s) -- construct a solution, starting at idx, of K nodes.
                ant'     = checkBest s ant bestAnt scoredata      -- keep it?

startNode :: Ant -> Edge -> Ant
-- give a new ant -- a place to start
startNode (a,b,c) e = (e, b, c)

makeNewAnt :: Edge -> Ant 
makeNewAnt e = (startNode newAnt e)

------------------------------------------------------------------------------------------------------------------

constructSoln :: Double -> Digraph -> (Ant,State) -> (Ant,State)
-- sample a sequence of nodes in the digraph.
constructSoln es dcg (ant,s) 
              | es < 1     = (ant,s) -- complete ant here ... add in the edge to go home.
              | es == 1    = constructSoln (es-1) dcg (getLastEdge s dcg ant)
              | otherwise  = constructSoln (es-1) dcg (getNextEdge s dcg ant)


getNextEdge :: State -> Digraph -> Ant -> (Ant,State)
getNextEdge s dcg ant = (ant', s{g=g'})
            where edgies        = possibleEdges ant dcg     :: Digraph    -- given where we are ... *
                  probabilities = edgeProb s edgies         :: [Double]   -- get probabilities for each possible edge 
                  cumulaProbs   = cumulate probabilities [] :: [Double]   -- make it a cumulative distribution
                  (dice, g')    = doubles 1 ([], (g s))     :: ([Double], StdGen) -- our random num
                  edgeIndex     = getIndex 0 ((head dice), cumulaProbs) :: Int    -- choose the edge
                  chosenEdge    = maybeSelect edgies edgeIndex :: Dedge   -- the chosen digraph edge
                  ant' = ((newNode ant chosenEdge), (chosenEdge : (solution ant)), 0.0)

-- final choice needs to lead the ant back to the starting node. --
getLastEdge :: State -> Digraph -> Ant -> (Ant,State)
getLastEdge s dcg ant = (ant', s{g=g'})
            where edgies        = possibleLastEdges ant dcg    :: Digraph   -- given where we are ... *
                  probabilities = edgeProb s edgies :: [Double]  -- get probabilities for each possible edge 
                  cumulaProbs   = cumulate probabilities []    :: [Double]  -- make it a cumulative distribution
                  (dice, g')    = doubles 1 ([], (g s))        :: ([Double], StdGen)  -- our random num
                  edgeIndex     = getIndex 0 ((head dice), cumulaProbs) :: Int   -- choose the edge
                  chosenEdge    = maybeSelect edgies edgeIndex :: Dedge
                  ant'          = ((newNode ant chosenEdge), (chosenEdge : (solution ant)), 0.0)

-- if K == 2, and the ant picks a self loop
-- then when K == 1, the ant can not find a path "back" to where it is .. (it's already traversed)
-- and the ant dies.

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- edgies = possibleEdges ant dcg

-- what edges can the ant see? -- Should be all of them!! -- That haven't already been walked on.
possibleEdges ::  Ant -> Digraph -> Digraph
possibleEdges ant dcg = filter (\x -> and [(connects ant x), (notEndingAtStart x), (notTraversed ant x)]) dcg
--possibleEdges ant dcg = filter (\x -> and [(notTraversed ant x), (connects ant x), (notEndingAtStart x)]) dcg 

-- what edges can the ant see? -- have to get back to the start.
-- ant either has been nowhere ..
-- or it has a path.
possibleLastEdges ::  Ant -> Digraph -> Digraph
possibleLastEdges (n1,es,s) dcg 
                  | es == []  = filter (\x -> loopsBack x) dcg -- no existing path
                  | otherwise = filter (\x -> and [(connects (n1,es,s) x),
                                                   (leadsBack n1 (last es) x),
                                                   (notTraversed (n1,es,s) x)]) dcg

-- has the ant traversed edge e?
-- this is where we might limit the memory of the ant..
-- type Ant = (Edge, Digraph, Double)  -- Current Location, Path, Score
notTraversed :: Ant -> Dedge -> Bool
notTraversed (ed, dig, sc) e = notInSoln e dig 

-- type Dedge = (Edge,Edge,Double,Double)    

notInSoln :: Dedge -> [Dedge] -> Bool
notInSoln e1 [] = True
notInSoln e1 (e2:es)
           | e1 == e2 = False 
           | otherwise = notInSoln e1 es

notEndingAtStart :: Dedge -> Bool -- Dedge (going-to)-(coming-from)
-- don't go back to the start until the end.  
notEndingAtStart (a1,b1,c1,d1) 
  | a1 == (("Start1","Start2",1.0,1.0):: Edge) = False
  | otherwise = True
                
-- does this ant sit on a node that connects to this edge?
-- dedges are going to be arranged:  (going-to)-(coming-from)
connects :: Ant -> Dedge -> Bool
connects (n1, es, s) (a,b,c,d)
         | n1 == b = True       -- why the c? -- && c > 0 = True
         | otherwise = False

--                       ax               bx
-- last is: ... (("r","s",0.99,0.5),("Start1","Start2",1.0,1.0)
--                       iy
-- want:        (("Start1","Start2",1.0,1.0),("i","j",0.99,0.5),0.0,0.5)
-- does this edge connect back to the given node?  --  (going-to)-(coming-from)
leadsBack :: Edge -> Dedge -> Dedge -> Bool   -- (i,j,_,_)        - x
leadsBack e1 (ax,bx,cx,dx) (iy,jy,ky,ly) -- (coming-from) - (going-to)
  | bx == iy && e1 == jy = True 
  | otherwise = False

-- do the two edges share a vertice?
shareVert :: Edge -> Edge -> Bool
shareVert (ax,bx,cx,dx) (ey,fy,gy,hy) 
          | ax == ey = True
          | ax == fy = True
          | bx == ey = True
          | bx == fy = True
          | otherwise = False

-- A self edge.
loopsBack :: Dedge -> Bool
loopsBack (ax,bx,cx,dx) 
          | ax == bx = True
          | otherwise = False

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- probabilities = List.map edgeProb edgies  
-- type Edge = (Node, Node, Double, Double)  -- From, To, EdgeType, EdgeWeight 
-- type Dedge = (Edge,Edge,Double,Double)    -- From, To, InfoCanPass, Pheromone
-- type Digraph = [Dedge]

-- list of probabilities of taking each edge
edgeProb :: State -> [Dedge] -> [Double]
edgeProb s es = map (\x -> x / sumscore) edgescore
         where edgescore = map (scoreFun s) es
               sumscore  = sum edgescore

-- the numerator of the probability function  
scoreFun :: State -> Dedge -> Double
scoreFun s ((n1,n2,d1,p1),(n3,n4,d2,p2),x,y) = ((d1*d2) ** (alph s)) * (y ** (beta s))

-- to produce a cumulative probabiltiy function
cumulate :: [Double] -> [Double] -> [Double]
cumulate []     ys = ys
cumulate (x:xs) [] = cumulate xs [x]
cumulate (x:xs) ys = (cumulate xs (ys ++ [(last ys) + x]))

-- which edge to randomly select
getIndex :: Int -> (Double, [Double]) -> Int
getIndex n (die, probs) 
  | probs == [] = (-1)
  | die == 0.0  = (-1)
  | die < (head probs) = n
  | otherwise          = getIndex (n+1) (die, (tail probs)) 

-- The ant might be stuck
maybeSelect x y 
  | y == -1 = deadDedge -- error, no possible edges --
  | otherwise = x !! y

newNode :: Ant -> Dedge -> Edge
newNode ant theEdge     
        | (node ant) == (firstist theEdge) = (secondist theEdge)
        | otherwise                        = (firstist theEdge)

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

-- Given the current solution, encoded as an ant, and the best ant, and the 
-- Directed Cyclic Graph (dcg) -- we score the ant ant return the better.
checkBest :: State -> Ant -> Ant -> ScoreData -> Ant
checkBest s currSolnAnt bestSolnAnt enm = x
          where trimmedAnt = trimAnt s currSolnAnt
                currSolnScore = score s trimmedAnt enm 
                currentAnt    = scoreTransfer currSolnScore currSolnAnt
                x             = bestAnt s currentAnt bestSolnAnt

trimAnt :: State -> Ant -> Ant
-- removes the start and end edges
trimAnt s (n1, es, scr) 
  | (k s) == 1 = (n1, (removeStartEdge es), scr)
  | otherwise = (n1, (init (tail es)), scr)
                
removeStartEdge :: Digraph -> Digraph
removeStartEdge es = newDigraph
                     where vertices   = map fstsnd es :: [[Edge]]
                           uVerts     = filter (\x -> noStart x) $ nub $ concat vertices ::  [Edge]   
                           newDigraph = [((head uVerts), (head uVerts), 1.0, 1.0)] :: Digraph


noStart :: Edge -> Bool
noStart e
  | e == ("Start1","Start2",1.0,1.0) = False
  | otherwise = True

scoreTransfer :: Ant -> Ant -> Ant
-- need the full path back
scoreTransfer from (toNode, toSoln, toScore) = (toNode, toSoln, (fitness from))

-- First time around, we need to get rid of the random ant.
bestAnt :: State -> Ant -> Ant -> Ant
bestAnt s currSolnAnt bestSolnAnt
        | (fitness bestSolnAnt) == 0.0 = currSolnAnt
        | (fitness currSolnAnt) > (fitness bestSolnAnt) = currSolnAnt
        | otherwise = bestSolnAnt

-------------------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------------------------------------
-- update solutions -- 
-------------------------------------------------------------------------------------------------------------------------
--type SolnSet = [(Double, [Edge])] -- best-ever, restart-best, iteration-best, 

updateSolns :: Ant -> State -> State
updateSolns ant s = s'
                where path   = solution ant :: Digraph
                      thisScore  = fitness ant  :: Double
                      s'     = checkUpdate s thisScore path

-- The head of solns is the best ever ... over all restarts
-- The middle is the best for this restart
-- The current best for this iteration of ants
checkUpdate :: State -> Double -> Digraph -> State
checkUpdate s thisScore path 
    -- new best so far
  | (fst (bestEver s)) <= thisScore && (fst (bestRest s) <= thisScore) = 
    s {bestEver=(thisScore, path), bestRest=(thisScore, path), bestIter=(thisScore, path)}
    -- best since the last restart
  | (fst (bestRest s) <= thisScore) = s {bestRest=(thisScore, path), bestIter=(thisScore, path)} 
    -- best this ant iteration
  | otherwise = s{bestIter = (thisScore, path)}

resetSol :: State -> State
-- reset just the restart and iteration bests
resetSol s = s {bestRest=(0.0, []), bestIter=(0.0, [])}

-------------------------------------------------------------------------------------------------------------------------
-- update Network --
-------------------------------------------------------------------------------------------------------------------------

updateDigraph :: Digraph -> State -> Digraph 
updateDigraph digraph s = updateDigraph' digraph s (proportional s) []

updateDigraph' :: Digraph -> State -> (Double,Double) -> Digraph -> Digraph
updateDigraph' [] s ab newgraph = newgraph
updateDigraph' (d:digraph) s ab newgraph 
  = updateDigraph' digraph s ab ((edgeupd s ab d) : newgraph )

-- take edge 
-- is the edge in the restart best or the iteration best or both or neither

edgeupd :: State -> (Double,Double) -> Dedge -> Dedge
edgeupd s (iterp,restp) (a1,b1,c1,d1) = (a1,b1,c1,zzz) -- add to the pheromone
                          where (inIter,inRest) = inSolutions s (a1,b1,c1,d1) 
                                deposited = inIter * iterp + inRest * restp
                                zzz = bounded $ d1 + (evap s) * (deposited - d1) 

proportional :: State -> (Double,Double)
-- (iteration best, restart best) --
proportional s
  | cval < 0.4                = (0.0,1.0)  -- nearing convergence, use restart best
  | cval >= 0.4 && cval < 0.6 = (0.3331, 0.6669)
  | cval >= 0.6 && cval < 0.8 = (0.6669, 0.3331)
  | cval >= 0.8               = (1.0, 0.0) -- just started out, use iteration best
  where cval = (c s)
                
inSolutions :: State -> Dedge -> (Double, Double)
-- edge is ... (in iteration best, in restart best)
inSolutions s dxy  
  | inRest && inIter = (1.0,1.0)
  | inRest           = (0.0,1.0)
  | inIter           = (1.0,0.0)
  | otherwise        = (0.0,0.0)
  where inRest = edgeElem dxy (snd (bestRest s)) -- needs to be in either direction
        inIter = edgeElem dxy (snd (bestIter s)) -- to update both edges

edgeElem :: Dedge -> [Dedge] -> Bool
-- doesn't matter what order the two edges are..
edgeElem d1 [] = False
edgeElem (a1,b1,c1,d1) ((e1,f1,g1,h1):ds)   
  | a1 == e1 && b1 == f1 = True
  | b1 == e1 && a1 == f1 = True
  | otherwise = edgeElem (a1,b1,c1,d1) ds

bounded :: Double -> Double
bounded b 
        | b < 0.001 = 0.001
        | b > 0.999 = 0.999
        | otherwise = b

-------------------------------------------------------------------------------------------------------------------------
-- check convergence --
-------------------------------------------------------------------------------------------------------------------------

-- will they be in the same order? -- 

convergence :: State -> Digraph -> State
-- all of the pheromone values are close to 0 or 1
convergence s digraph = s {c = convergeFactor, ptoKratio = ((sum doublelist) / (2*(0.999-0.001)*((k s)+1)))}
            where normfactor = (0.999-0.001) * (fromIntegral (length digraph) :: Double)
                  doublelist = pheroVals digraph  
                  convergeFactor = 1 - (2 * (((sum (map convNum doublelist)) / normfactor)  - 0.5))
                  
--trace (show (take 10 (pheroVals digraph)))

convNum :: Double -> Double
convNum x1 = maximum [0.999-x1, x1-0.001]

pheroVals :: Digraph -> [Double]
pheroVals digraph = map (\(a1,b1,c1,d1) -> d1) digraph

-------------------------------------------------------------------------------------------------------------------------
-- Local Optimization --
-------------------------------------------------------------------------------------------------------------------------

localOptimize :: State -> Digraph -> ScoreData -> Ant -> Ant
localOptimize s digraph scoredata ant
              | (local s) == 1 = localOptimize' s digraph scoredata ant
              | otherwise = ant

localOptimize' :: State -> Digraph -> ScoreData -> Ant -> Ant
-- we want alternate solutions
localOptimize' s digraph scoredata (node, thissoln, thisscore) = theNewAnt
              where newSolns  = makeSolns thissoln digraph :: [Digraph]
                    newAnts   = map (\x -> (node, x, 0.0)) newSolns :: [Ant]
                    newScores = map (\x -> (fitness (score s x scoredata))) newAnts :: [Double]
                    (soln,idx) = maybeUpdate newSolns newScores thissoln thisscore
                    theNewAnt = (node, soln, idx) 

maybeUpdate :: [Digraph] -> [Double] -> Digraph -> Double -> (Digraph,Double)
maybeUpdate [] [] a b = (a,b)  -- just keep the ant
maybeUpdate (asoln:solns) (ascore:scores) oneSoln oneScore 
           | ascore > oneScore =  (asoln, ascore)  -- this one's better
           | otherwise = maybeUpdate solns scores oneSoln oneScore

--trace ("Local Opt. " ++ (show ascore) ++ "   " ++ (show oneScore)) 

-------------------------------------------------------------------------------------------------------------------------
-- END OF LINE --
-------------------------------------------------------------------------------------------------------------------------
-- ************************************************************************************************************************

