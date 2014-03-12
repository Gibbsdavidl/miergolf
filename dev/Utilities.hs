-------------------------------------------------------------------------------------------------------------------------
-- Influence Maximization Problem using hypercube-minmax ant optimization ----------------------------------------------
-- for problems on biological pathway graphs 
-------------------------------------------------------------------------------------------------------------------------

-- Input:
-- directed pathways with edge weights --

-- Output:
-- approximate best subset of size K, with the greatest influence on graph, by diffusion models.
-------------------------------------------------------------------------------------------------------------------------


-- {-# OPTIONS -Wall #-} 

module Utilities (EdgeList,
                  Node,
                  Edge,
                  Network,
                  Idx,
                  Digraph,
                  Ant,
                  doubles,
                  emptyEdge,
                  newAnt,
                  deadDedge,
                  first,
                  second,
                  third,
                  firstist,
                  secondist,
                  thirdist,
                  lastist,
                  fstsnd,
                  node,
                  solution,
                  fitness,
                  listNodes,
                  canPass,
                  digraphSize,
                  emptySM,
                  emptyDigraph,
                  emptyIdx,
                  emptySV
                  ) where
----------------------------------------------------------------------------------------------------
-- The Network data type and functions on the network ----------------------------------------------
----------------------------------------------------------------------------------------------------

import System.Random
import Numeric.LinearAlgebra
import qualified Data.IntMap.Strict as IM (IntMap, lookup, map, unionWith, empty, size,
                                           intersectionWith, fromList, findWithDefault, foldl')
import qualified Data.Map as Map
import Sparse

type EdgeList = [[String]]

type Node = String

--   Edge =  Node1 Node2 Weight  Extra
type Edge = (Node, Node, Double, Double)  -- From, To, EdgeType, EdgeWeight 
type Network = [Edge]

type Path = [Int] -- a set of indices into the sparse matrix.

type Idx = Map.Map Edge Int
  
type Digraph = (Idx, SM Double)

type Ant = (Int, Path, Double)  -- Current Location, Path, Score



-------------------------------------------------------------------------------------------------------------------------
-- utility functions --
-------------------------------------------------------------------------------------------------------------------------

-- generate lists of random doubles
doubles :: Int -> ([Double], StdGen) -> ([Double], StdGen)
doubles 0 (xs, g) = (xs, g)
doubles n (xs, g) = doubles (n-1) ((z:xs), g')
        where (z, g') = randomR (0.0, 1.0) g 
              
double1 :: Int -> StdGen -> [Double]
double1 n g = fst $ doubles n ([],g)

emptyEdge = (-1)
newAnt = (emptyEdge, [], 0.0) :: Ant
deadDedge = 0 :: Int
emptyIdx = Map.fromList [(("xa","xb",0.0,0.0),0)] :: Map.Map Edge Int
emptySV  = SV (IM.empty)
emptySM  = SM (0, IM.fromList [   (0,SV (IM.fromList [(0,(0.0,0.0))])) ])
emptyDigraph = (emptyIdx, emptySM)

digraphSize (x, (SM (_, m))) = show $ IM.size m

first :: (t, t1, t2) -> t
first (a,_,_) = a
second :: (t, t1, t2) -> t1
second (_,b,_) = b
third :: (t, t1, t2) -> t2
third (_,_,c) = c

firstist :: (t,t1,t2,t3) -> t
firstist (a,_,_,_) = a

secondist :: (t,t1,t2,t3) -> t1
secondist (_,a,_,_) = a

thirdist :: (t,t1,t2,t3) -> t2
thirdist (_,_,a,_) = a

lastist :: (t,t1,t2,t3) -> t3
lastist (_,_,_,a) = a

fstsnd :: (t2, t2, t, t1) -> [t2]
fstsnd (a,b,_,_) = [a,b] 

node :: Ant -> Int
node (a,_,_) = a

solution :: Ant -> Path
solution (_,a,_) = a

fitness :: Ant -> Double
fitness (_,_,c) = c

-- returning a list of the nodes within a network
listNodes :: Network -> [String]
listNodes l = listNodes' l []

listNodes' :: Network -> [String] -> [String]
listNodes' [] ll = ll
listNodes'  l  ll = listNodes' (tail l) ( (fstsnd (head l)) ++ ll)

canPass :: (Edge, Edge, Double, Double) -> Bool
canPass (_,_,a,_)
  | a > 0.0  = True
  | otherwise = False


-------------------------------------------------------------------------------------------------------------------------
-- END OF LINE --
-------------------------------------------------------------------------------------------------------------------------

-- ************************************************************************************************************************
