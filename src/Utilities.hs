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
                  Dedge,
                  Digraph,
                  Ant,
                  ScoreData,
                  doubles,
                  emptyEdge,
                  newAnt,
                  deadDedge,
                  first,
                  second,
                  third,
                  firstist,
                  secondist,
                  lastist,
                  fstsnd,
                  node,
                  solution,
                  fitness,
                  listNodes,
                  canPass,
                  ) where
----------------------------------------------------------------------------------------------------
-- The Network data type and functions on the network ----------------------------------------------
----------------------------------------------------------------------------------------------------

import System.Random
import Numeric.LinearAlgebra

type EdgeList = [[String]]

type Node = String

--   Edge =  Node1 Node2 Weight  Extra
type Edge = (Node, Node, Double, Double)  -- From, To, EdgeType, EdgeWeight 
type Network = [Edge]

type Dedge = (Edge,Edge,Double,Double) -- From, To, InfoCanPass, Pheromone
type Digraph = [Dedge]

type Ant = (Edge, Digraph, Double)  -- Current Location, Path, Score

type ScoreData = (Network, Matrix Double)



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

emptyEdge = ("%x%","%y%",0.0,0.0) :: Edge
newAnt = (emptyEdge, [], 0.0) :: Ant
deadDedge = (("%x&","%y&",0.0,0.0),("&z%","&w%",0.0,0.0),0.0,0.0) :: Dedge

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

lastist :: (t,t1,t2,t3) -> t3
lastist (_,_,_,a) = a

fstsnd :: (t2, t2, t, t1) -> [t2]
fstsnd (a,b,_,_) = [a,b] 

node :: Ant -> Edge
node (a,_,_) = a

solution :: Ant -> Digraph
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
