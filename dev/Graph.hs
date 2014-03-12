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

-- Idx :: intmap of (Ints, Edges)

-- let a = SM (2, 
--            fromList [
--                     (0, SV (fromList [(0, 4), (1, 1)])),  -- each row is a list of tuples  
--                     (1, SV (fromList [(0, 1), (1, 3)]))   -- (zipped with index) 
--                     ]
--            ) :: SM Double



module Graph
(initGraphData
,initGraphDataT
,constructDigraph
,doesPass
,resetNet
,initSMList
--,resetDigraph
,showSize) where 

import Sparse
import ProgramState (State,emptyState,graphfile)
import Utilities
import Numeric.LinearAlgebra
import Data.Maybe
import System.IO.Unsafe
import qualified Data.IntMap.Strict as IM
import qualified Data.Map.Strict as Map

initGraphDataT :: Maybe State -> String -- (State, Digraph)
initGraphDataT s 
  | isNothing s = "nothing"
  | isJust s    = something ++ "\n" ++ (prettyDigraphPrint digraph)
                  where digraph  = constructDigraph (fromJust s)
                        something = digraphSize digraph
                  
initGraphData :: Maybe State -> (State, Digraph)
initGraphData s 
  | isNothing s = (emptyState, emptyDigraph)
  | isJust s    = ((fromJust s), digraph)
                  where digraph  = constructDigraph (fromJust s)

constructDigraph :: State -> Digraph
constructDigraph s = buildDigraph $ buildNetwork $ Prelude.map words (lines (unsafePerformIO (readFile (graphfile s))))

-- build the edge index --
buildNetwork :: EdgeList -> Idx
buildNetwork l = buildNetwork' 1 l (Map.fromList [(("Start1", "Start2", 1.0, 1.0), 0)])

buildNetwork' :: Int -> EdgeList -> Idx -> Idx
buildNetwork' i l n
  | length l == 0 = n -- the start node -- init the IntMap here
  | otherwise = buildNetwork' (i+1) (tail l) (Map.insert ((hl !! 0), (hl !! 1), (read (hl !! 2) :: Double), 0.5) i n)
    where hl  = head l

-- build the sparse matrix using the Idx--
buildDigraph :: Idx -> Digraph
buildDigraph net = zymorph
  where n       = sideSize net     -- get the size of one side of the matrix
        smlist  = initSMList n     -- our empty sparse matrix
        epairs  = allEdgePairs net -- for each potential entry in the matrix
        smlist' = addEntry epairs net smlist  -- if infoPasses, then it's added to the matrix
        smlist''= filter (\(a, m) -> emptySV m) smlist'
        zymorph = (net, SM (n, IM.fromList smlist'')) 
                 
-- Idx :: intmap of (Ints, Edges)

-- let a = SM (2, 
--            fromList [
--                     (0, SV (fromList [(0, 4), (1, 1)])),  -- each row is a list of tuples  
--                     (1, SV (fromList [(0, 1), (1, 3)]))   -- (zipped with index) 
--                     ]
--            ) :: SM Double

sideSize :: Idx -> Int
-- The size of one side of the matrix --
sideSize ime = length $ listNodes $ (Map.keys ime) -- get the nodes out of the list of edges


initSMList :: Int -> [(Int, SV Double)]
-- build a list, a SV for each node in the Idx list ... these should be filtered later
initSMList n = zip [0 .. n] (replicate n anEmptySV)

-- let x = initSMList 10
-- let y = SM (10, IM.fromList x)

allEdgePairs :: Idx -> [(Edge,Edge)]
allEdgePairs ime = [(x,y) | x <- (Map.keys ime), y <- (Map.keys ime)]

addEntry :: [(Edge,Edge)] -> Idx -> [(Int, SV Double)] -> [(Int, SV Double)]
-- take the pair, if info passes, add the wt*wt and pheromone
-- to list position given by e1 and posiiton
addEntry [] idx smlist = smlist
addEntry ((e1,e2):es) idx smlist 
  | doesPass e1 e2 = addEntry es idx (insertD e1 e2 idx smlist)
  | otherwise = addEntry es idx smlist

insertD :: Edge -> Edge -> Idx -> [(Int, SV Double)] -> [(Int, SV Double)]
insertD e1 e2 idx smlist = smlist'
  where i1 = fromJust $ Map.lookup e1 idx
        i2 = fromJust $ Map.lookup e2 idx
        wt = (thirdist e1) * (thirdist e2)
        smlist' = Prelude.map (listInsert i1 i2 wt) smlist
          
listInsert :: Int -> Int -> Double -> (Int, SV Double) -> (Int, SV Double)
listInsert i1 i2 wt (n, SV m)
  | i1 == n = (n, SV (IM.insert i2 (wt, 0.5) m))
  | otherwise = (n, SV m)

-- for building the matrix :: (coming from) - (going to)
doesPass :: Edge -> Edge -> Bool -- can info pass from this edge to that edge?
doesPass (_,a,_,_) (b,_,_,_)
  | a == b    = True
  | otherwise = False

-- type Edge = (Node, Node, Double, Double)  --

resetNet :: Fractional t4 => [(t1, t2, t3, t)] -> [(t1, t2, t3, t4)]
resetNet network = resetNet' network []

resetNet' :: Fractional t4 => [(t1, t2, t3, t)] -> [(t1, t2, t3, t4)] -> [(t1, t2, t3, t4)]
resetNet' [] x = x
resetNet' ((a,b,c,_):net) x = resetNet' net ( (a,b,c,0.5):x )

-- type Dedge = (Edge,Edge,Double,Double) -- From, To, InfoCanPass, Pheromone

--resetDigraph :: Digraph -> Digraph
--resetDigraph digraph = resetNet digraph

subsetSM :: [Int] -> SM a -> [SV a]
subsetSM xs (SM (_, m)) = catMaybes $ Prelude.map (\x -> (x `IM.lookup` m)) xs

subsetSV :: [Int] -> SV a -> SV a
subsetSV [] (SV v) = (SV v)
subsetSV (x:xs) (SV v) = subsetSV xs (SV (IM.delete x v))

showSize :: SM a -> String
showSize (SM (_, m)) = show $ IM.size m



-------------------------------------------------------------------------------------------------------------------------
-- END OF LINE --
-------------------------------------------------------------------------------------------------------------------------

-- ************************************************************************************************************************

