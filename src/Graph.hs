-------------------------------------------------------------------------------------------------------------------------
-- Influence Maximization Problem using hypercube-minmax ant optimization ----------------------------------------------
-- for problems on biological pathway graphs 
-------------------------------------------------------------------------------------------------------------------------


-- main = getArgs >>=  printResult . antopt . buildGraphs . buildState
--                                               ***                                                  
-- Input:
-- State records

-- Output:
-- Tuple, holding the state records, and the graphs.
-------------------------------------------------------------------------------------------------------------------------

module Graph
(initGraphData
,constructDigraph
,doesPass
,resetNet
,resetDigraph) where 

import Utilities
import Numeric.LinearAlgebra
import ProgramState
import Data.Maybe
import InfoFlow
import System.IO.Unsafe


initGraphData :: Maybe State -> (State, Digraph, ScoreData)
initGraphData s 
  | isNothing s = (emptyState, [], ([], ident 2))
  | isJust s    = ((fromJust s), digraph, scoredata)
                  where digraph   = constructDigraph (fromJust s)
                        scoredata = transformer (fromJust s) digraph

constructDigraph :: State -> Digraph
constructDigraph s = buildDigraph $ buildNetwork $ map words (lines (unsafePerformIO (readFile (graphfile s))))

buildNetwork :: EdgeList -> Network
buildNetwork l = buildNetwork' l []

buildNetwork' :: EdgeList -> Network -> Network
buildNetwork' l n 
  | length l == 0 = ("Start1", "Start2", 1, 1) : n -- the start node.
  | otherwise = buildNetwork' (tail l) (((hl !! 0), (hl !! 1), (read (hl !! 2) :: Double), 0.5) : n)
    where hl  = head l

buildDigraph :: Network -> Digraph
buildDigraph net = [(a,b,(doesPass a b),0.5) | a <- net, b <- net] -- we get dedges

-- for building paths dedges go :: (going-to) - (coming-from)
doesPass :: Edge -> Edge -> Double -- can info pass from this edge to that edge?
doesPass (_,a,_,_) (b,_,_,_)
  | a == b = 1.0
  | otherwise = 0.0

-- type Edge = (Node, Node, Double, Double)  --

resetNet :: Fractional t4 => [(t1, t2, t3, t)] -> [(t1, t2, t3, t4)]
resetNet network = resetNet' network []

resetNet' :: Fractional t4 => [(t1, t2, t3, t)] -> [(t1, t2, t3, t4)] -> [(t1, t2, t3, t4)]
resetNet' [] x = x
resetNet' ((a,b,c,_):net) x = resetNet' net ( (a,b,c,0.5):x )

-- type Dedge = (Edge,Edge,Double,Double) -- From, To, InfoCanPass, Pheromone

resetDigraph :: Digraph -> Digraph
resetDigraph digraph = resetNet digraph


-------------------------------------------------------------------------------------------------------------------------
-- END OF LINE --
-------------------------------------------------------------------------------------------------------------------------

-- ************************************************************************************************************************

