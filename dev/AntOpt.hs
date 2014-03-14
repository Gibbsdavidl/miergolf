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
,edgeProb
,cumulate
,getIndex
,generateSoln
--,localOptimize
,convergence
,updateSolns
,proportional
,inSolutions
,bounded
) where


import qualified Data.IntMap.Strict as IM
import Data.Maybe as Maybe (fromJust)
import ProgramState
import Utilities
import Graph
import Sparse
-- import LocOpt
import System.Random
import Debug.Trace
import Numeric.LinearAlgebra
import System.IO.Unsafe
import Data.List (nub)

optimize :: (State, Digraph) -> String
optimize (state,(idx, sparse)) 
  | isEmptyState state = "ERROR"  -- stop now!
  | otherwise = showState (minmax state idx sparse)  -- OK, let's compute
 
------------------------------------------------------------------------------------------
-- The particular optimization routine -- MinMax Ant Optimization
-----------------------------------------------------------------------------------------
minmax :: State -> Idx -> SM Double -> State
minmax s idx sparse
  | doneState s = s  -- all done? , return program state
  | converged s = minmax (newRun s) (resetIdx idx) sparse 
  | otherwise   = minmax s' idx' sparse -- or do another iteration --
    where (s', idx') = runColony s idx sparse

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
      
runColony :: State -> Idx -> SM Double -> (State, Idx)   -- changed states: StdGen,SolnSet,Runs
-- An Ant iteration: run colony needs to return a digraph and a state.
runColony s idx sparse = (s3, idx')
  where (ant, s')  = generateSoln s idx sparse -- generate solutions -- keep the best
        -- ant2       = localOptimize s' digraph scoredata ant
        s2         = updateSolns ant s'       -- is it better than the previous soln?
        idx'       = updateIdx idx s2          -- update the pheromone levels
        s3         = convergence s2 idx'   -- check for convergence


------------------------------------------------------------------------------------------------
-- compute solutions -- 
-------------------------------------------------------------------------------------------------------------------------

generateSoln :: State -> Idx -> SM Double -> (Ant, State)
-- For i = 1 to m ... create an ant ... construct a solution.
--     if the solution is better then a held solution ... keep it.
generateSoln s idx sparse = (genSoln (ants s) s idx probs sparse newAnt)
  where ks = ([0..((IM.size idx) - 1)])
        probs = edgeProb s ks idx :: [Double]

genSoln :: Int -> State -> Idx -> [Double] -> SM Double -> Ant -> (Ant,State)
-- generate a solution for each ant --
genSoln 0    s idx ps sparse bestAnt = (bestAnt, s)
genSoln curr s idx ps sparse bestAnt = genSoln (curr-1) s' idx ps sparse ant'
          where newant   = makeNewAnt (-1)                       
                (ant,s') = constructSoln (k s) idx (newant, s, (IM.keys idx), ps)
                ant'     = checkBest s' ant bestAnt sparse idx      -- keep it?

makeNewAnt :: Int -> Ant 
makeNewAnt e = (startNode newAnt e)

startNode :: Ant -> Int -> Ant
-- give a new ant -- a place to start
startNode (a,b,c) e = (e, b, c)


------------------------------------------------------------------------------------------------------------------

constructSoln :: Double -> Idx -> (Ant,State,[Int],[Double]) -> (Ant,State)
-- sample a set of nodes in the digraph.
constructSoln es idx (ant,s,ks,ps) 
              | es < 1     = trace (show ant) (ant,s) -- complete ant here ... add in the edge to go home.
              | otherwise  = constructSoln (es-1) idx (getNext s idx ks ps ant) -- ks, ps will change...


getNext :: State -> Idx -> [Int] -> [Double] -> Ant -> (Ant,State, [Int],[Double])
-- ks are the set of items left to choose from
getNext s idx ks ps ant = (ant', s{g=g'}, ks', ps')
        -- Maybe first filter out ks that are not interesting .. --
  where cumulaProbs   = cumulate ps [] :: [Double]                         -- cumulative distribution
        (dice, g')    = doubles 1 ([], (g s))     :: ([Double], StdGen)    -- our random num
        i             = getIndex 0 ((head dice), cumulaProbs) :: Int -- choose the edge
        k             = ks !! i
        (ks',ps')     = removeChosenK 0 i ks ps ([], [])
        ant'          = (k, (k : (solution ant)), 0.0)

removeChosenK :: Int -> Int -> [Int] -> [Double] -> ([Int],[Double]) -> ([Int],[Double])
removeChosenK i j (ki:ks) (pi:ps) (x,y)  
  -- get rid of the first one
  | j == 0 = (ks, ps)
  -- keep some, get rid of ki and pi, but keep rest 
  | i == j = (x ++ ks, y ++ ps)
  | otherwise = removeChosenK (i+1) j ks ps (x++[ki], y++[pi])

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

-- list of probabilities of taking each edge
-- We could probably make this more efficient since many of the entries will ahve the sampe
-- probabiliies --- if they are not being selected..
edgeProb :: State -> [Int] -> Idx -> [Double]
edgeProb s ks idx = map (\x -> x / sumscore) edgescore                           -- ASSUMING --
         where edgescore = [(compScore s i idx) | i <- ks] :: [Double] -- SHOULD BE IN ORDER OF KS --
               sumscore  = foldl (+) 0 edgescore

compScore :: State -> Int -> Idx -> Double
compScore s i idx = scoreFun s (Maybe.fromJust (IM.lookup i idx))

-- the numerator of the probability function  
scoreFun :: State -> Edge -> Double
scoreFun s (n1,n2,w1,p1,x) = (w1 ** (alph s)) * (p1 ** (beta s))

-- to produce a cumulative probabiltiy function -- retain ordering!
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

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

-- Given the current solution, encoded as an ant, and the best ant, and the 
-- Directed Cyclic Graph (dcg) -- we score the ant ant return the better.
checkBest :: State -> Ant -> Ant -> SM Double -> Idx -> Ant
checkBest s currSolnAnt bestSolnAnt sparse idx = trace ("x: " ++ show x) x
          where currentAnt = score s currSolnAnt sparse idx
                x          = bestAnt s currentAnt bestSolnAnt

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
                where path      = solution ant :: [Int]
                      thisScore = fitness  ant :: Double
                      s'        = checkUpdate s thisScore path

-- The head of solns is the best ever ... over all restarts
-- The middle is the best for this restart
-- The current best for this iteration of ants
checkUpdate :: State -> Double -> [Int] -> State
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

updateIdx :: Idx -> State -> Idx
updateIdx idx s = updateIdx' idx s (proportional s)

updateIdx' :: Idx -> State -> (Double,Double) -> Idx
-- the best results are carried in the state
updateIdx' idx s (iterp,restp) = IM.fromList [(newIdx i idx s iterp restp) | i <- IM.keys idx] 
                                        
newIdx :: Int -> Idx -> State -> Double -> Double -> (Int, Edge)
-- is this edge in the solutions, iter or restart bests?
newIdx i idx s iterp restp = (i, newEdge)
  where (n1,n2,wt,ph,pi) = fromJust $ IM.lookup i idx
        (inIter, inRest) = inSolutions s i
        deposited = inIter * iterp + inRest * restp
        ph' = bounded $ ph + (evap s) * (deposited - ph)
        newEdge = (n1,n2,wt,ph',(inIter+inRest+pi)) :: Edge

proportional :: State -> (Double,Double)
-- (iteration best, restart best) --
proportional s
  | cval < 0.4                = (0.0,1.0)  -- nearing convergence, use restart best
  | cval >= 0.4 && cval < 0.6 = (0.3331, 0.6669)
  | cval >= 0.6 && cval < 0.8 = (0.6669, 0.3331)
  | cval >= 0.8               = (1.0, 0.0) -- just started out, use iteration best
  where cval = (c s)
                
inSolutions :: State -> Int -> (Double, Double)
-- edge is ... (in iteration best, in restart best)
inSolutions s i
  | inRest && inIter = (1.0,1.0)
  | inRest           = (0.0,1.0)
  | inIter           = (1.0,0.0)
  | otherwise        = (0.0,0.0)
  where inRest = elem i (snd (bestRest s)) -- needs to be in either direction
        inIter = elem i (snd (bestIter s)) -- to update both edges
        
bounded :: Double -> Double
bounded b 
        | b < 0.001 = 0.001
        | b > 0.999 = 0.999
        | otherwise = b

-------------------------------------------------------------------------------------------------------------------------
-- check convergence --
-------------------------------------------------------------------------------------------------------------------------

-- will they be in the same order? -- 

convergence :: State -> Idx -> State
-- all of the pheromone values are close to 0 or 1
convergence s idx = s {c = convergeFactor, ptoKratio = ((sum doublelist) / (2*(0.999-0.001)*((k s)+1)))}
            where normfactor = (0.999-0.001) * (fromIntegral (IM.size idx) :: Double)
                  doublelist = pheroVals idx  
                  convergeFactor = 1 - (2 * (((sum (map convNum doublelist)) / normfactor)  - 0.5))
                  
--trace (show (take 10 (pheroVals digraph)))

convNum :: Double -> Double
convNum x1 = maximum [0.999-x1, x1-0.001]

pheroVals :: Idx -> [Double]
pheroVals idx = Prelude.map (\(a1,b1,c1,d1,e1) -> d1) (IM.elems idx)

-------------------------------------------------------------------------------------------------------------------------
-- Local Optimization --
-------------------------------------------------------------------------------------------------------------------------
{-
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
-}
-------------------------------------------------------------------------------------------------------------------------
-- END OF LINE --
-------------------------------------------------------------------------------------------------------------------------
-- ************************************************************************************************************************

