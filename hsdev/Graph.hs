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
,constructDigraph
,doesPass
,resetIdx
,initSMList
,Digraph
,showSize) where 

import Sparse
import ProgramState (State,emptyState,graphfile)
import Utilities
import Numeric.LinearAlgebra
import Data.Maybe
import System.IO.Unsafe
import qualified Data.IntMap.Strict as IM

type Digraph = (Idx, SM Double)

emptyDigraph = (emptyIdx, emptySM)
prettyDigraphPrint :: Digraph -> String
prettyDigraphPrint (x,(SM (_, m))) = IM.showTree m

initGraphData :: Maybe State -> (State, Digraph)
initGraphData s 
  | isNothing s = (emptyState, emptyDigraph)
  | isJust s    = ((fromJust s), digraph)
                  where digraph  = constructDigraph (fromJust s)

constructDigraph :: State -> Digraph
constructDigraph s = buildDigraph $ buildNetwork $ Prelude.map words (lines (unsafePerformIO (readFile (graphfile s))))

-- build the edge index --
buildNetwork :: EdgeList -> Idx
buildNetwork l = buildNetwork' 0 l (IM.empty)

buildNetwork' :: Int -> EdgeList -> Idx -> Idx
buildNetwork' i l n
  | length l == 0 = n -- the start node -- init the IntMap here
  | otherwise = buildNetwork' (i+1) (tail l) (IM.insert i ((hl !! 0), (hl !! 1), (read (hl !! 2) :: Double), 0.5, 0) n)
    where hl  = head l

-- build the sparse matrix using the Idx--
buildDigraph :: Idx -> Digraph
buildDigraph idx = zymorph
  where n       = sideSize idx     -- get the size of one side of the matrix
        smlist  = initSMList n     -- our empty sparse matrix
        epairs  = allEdgePairs idx -- for each potential entry in the matrix
        smlist' = addEntry epairs idx smlist  -- if infoPasses, then it's added to the matrix
        smlist''= filter (\(a, m) -> emptySV m) smlist'
        zymorph = (idx, SM (n, IM.fromList smlist'')) 
                 
-- Idx :: intmap of (Ints, Edges)

-- let a = SM (2, 
--            fromList [
--                     (0, SV (fromList [(0, 4), (1, 1)])),  -- each row is a list of tuples  
--                     (1, SV (fromList [(0, 1), (1, 3)]))   -- (zipped with index) 
--                     ]
--            ) :: SM Double

sideSize :: Idx -> Int
-- The size of one side of the matrix --
sideSize ime = length $ (IM.keys ime) -- get the nodes out of the list of edges


initSMList :: Int -> [(Int, SV Double)]
-- build a list, a SV for each node in the Idx list ... these should be filtered later
initSMList n = zip [0 .. n] (replicate n anEmptySV)

-- let x = initSMList 10
-- let y = SM (10, IM.fromList x)

allEdgePairs :: Idx -> [(Int,Int)]
allEdgePairs ime = [(x,y) | x <- (IM.keys ime), y <- (IM.keys ime)] -- (Int, Int)

addEntry :: [(Int,Int)] -> Idx -> [(Int, SV Double)] -> [(Int, SV Double)]
-- take the pair, if info passes, add the wt*wt and pheromone
-- to list position given by e1 and posiiton
addEntry [] idx smlist = smlist
addEntry ((i1,i2):is) idx smlist 
  | checkPass i1 i2 idx = addEntry is idx (insertD i1 i2 idx smlist)
  | otherwise = addEntry is idx smlist

insertD :: Int -> Int -> Idx -> [(Int, SV Double)] -> [(Int, SV Double)]
insertD i1 i2 idx smlist = smlist'
  where e1 = fromJust $ IM.lookup i1 idx
        e2 = fromJust $ IM.lookup i2 idx
        wt = (thirdist e1) * (thirdist e2)
        smlist' = Prelude.map (listInsert i1 i2 wt) smlist
          
listInsert :: Int -> Int -> Double -> (Int, SV Double) -> (Int, SV Double)
listInsert i1 i2 wt (n, SV m)
  | i1 == n = (n, SV (IM.insert i2 wt m))
  | otherwise = (n, SV m)


-- for building the matrix :: (coming from) - (going to)
checkPass :: Int -> Int -> Idx -> Bool
checkPass i1 i2 idx = doesPass (fromJust (IM.lookup i1 idx)) (fromJust (IM.lookup i2 idx))

-- for building the matrix :: (coming from) - (going to)
doesPass :: Edge -> Edge -> Bool -- can info pass from this edge to that edge?
doesPass (_,a,_,_,_) (b,_,_,_,_)
  | a == b    = True
  | otherwise = False

-- type Edge = (Node, Node, Double, Double)  --

resetIdx :: Idx -> Idx
resetIdx idx = IM.fromList [(i, (resetEdge (fromJust (IM.lookup i idx)))) | i <- IM.keys idx]

resetEdge :: Edge -> Edge
resetEdge (ax,bx,cx,dx,ex) = (ax,bx,cx,0.5,ex)

subsetSM :: [Int] -> SM a -> [SV a]
subsetSM xs (SM (_, m)) = catMaybes $ Prelude.map (\x -> (x `IM.lookup` m)) xs

subsetSV :: [Int] -> SV a -> SV a
subsetSV [] (SV v) = (SV v)
subsetSV (x:xs) (SV v) = subsetSV xs (SV (IM.delete x v))

showSize :: SM a -> String
showSize (SM (_, m)) = show $ IM.size m

anEmptySV  = SV (IM.empty)
emptyIdx = IM.fromList [(0,("xa","xb",0.0,0.0,0.0))] :: IM.IntMap Edge
emptySM  = SM (0, IM.fromList [   (0,SV (IM.fromList [(0,(0.0))])) ])
digraphSize (x, (SM (_, m))) = show $ IM.size m


-------------------------------------------------------------------------------------------------------------------------
-- END OF LINE --
-------------------------------------------------------------------------------------------------------------------------

-- ************************************************************************************************************************

