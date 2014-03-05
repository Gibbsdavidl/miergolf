-------------------------------------------------------------------------------------------------------------------------
-- Influence Maximization Problem using hypercube-minmax ant optimization ----------------------------------------------
-- for problems on biological pathway graphs 
-------------------------------------------------------------------------------------------------------------------------

-- main = getArgs >>=  printResult . antOpt . buildGraphs . buildState
--                                                ***
-- Input:
-- To build the matrix ... need the digraph
-- To score ant solutions, need the ant containing the solution (a digraph).

-- Output:
-- A score.
-------------------------------------------------------------------------------------------------------------------------

module InfoFlow 
       (ScoreData
       ,transformer
       ,normItUpG
       ,fillInList
       ,updateList
       ,checkDiag
       ,replaceIJ
       ,replaceJ
       ,extractCols
       ,nm
       ,accumIndex
       ,makePTSandPST
       ,makePTT
       ,scoreGraph
       ,getScoreGraphMatrix
       ,score
       ,hasDeadEdge
       ) where

import Utilities
import ProgramState (State,damp,rx,tx)
import Numeric.LinearAlgebra
import Data.List (nub, elemIndex, (\\))
import Debug.Trace

transformer :: State -> Digraph -> ScoreData
transformer s dcg    = (uVerts, matrix')
  where vertices   = map fstsnd dcg                       :: [[Edge]]
        uVerts     = nub $ concat vertices                ::  [Edge]   -- list of edges!!!
        n          = length uVerts                        ::   Int
        alist      = replicate n $ take n (repeat 0.0)    :: [[Double]]
        filledlist =  fillInList uVerts dcg alist         :: [[Double]]
        normlist   = map (\x -> normItUpG s x) filledlist :: [[Double]]
        matrix'    = fromLists normlist :: Matrix Double

normItUpG :: State -> [Double] -> [Double]
normItUpG s doubleList
          | (sum doubleList) <= 0.0 = take (length doubleList) (repeat 0)
          | otherwise = map (\y -> (damp s)*(y/(sum doubleList))) doubleList 

fillInList :: [Edge] -> Digraph -> [[Double]] -> [[Double]]
fillInList edges [] mat          = mat
fillInList edges (dedge:dcg) mat = fillInList edges dcg (updateList edges mat dedge)

updateList :: [Edge] -> [[Double]] -> Dedge -> [[Double]]
updateList edges aList (a,b,c,d) = newList
  where i = maybe (-1) id $ elemIndex a edges -- where the a is in edge list
        j = maybe (-1) id $ elemIndex b edges -- where b is in the edge list
        v = checkDiag i j ((lastist a)*(lastist b))*c -- can info flow??
        newList = replaceIJ i j v aList
        
checkDiag :: Int -> Int -> Double -> Double
checkDiag i j d 
  | i == j = 0.0
  | otherwise = d
        
replaceIJ :: Int -> Int -> Double -> [[Double]] -> [[Double]]
replaceIJ i j newVal (x:xs)
     | i == 0 = (replaceJ j newVal x) : xs
     | otherwise = x:replaceIJ (i-1) j newVal xs
        
replaceJ n newVal (x:xs)
     | n == 0 = newVal:xs
     | otherwise = x:replaceJ (n-1) newVal xs

extractCols :: Element t => [Int] -> Matrix t -> Matrix t
extractCols l m = fromColumns $ extract (toColumns m) l
  where extract l' is = [l'!! i | i <- is]

nm :: Maybe Int -> Int
nm x = maybe (-1) id x

accumIndex :: [Int] -> [Int]
accumIndex xs = accumIndex' 0 xs []

accumIndex' :: Int -> [Int] -> [Int] -> [Int]
accumIndex' i [] ys = ys
accumIndex' i (x:xs) ys 
            | x > (-1)  = accumIndex' (i+1) xs (i:ys)
            | otherwise = accumIndex' (i+1) xs ys

-- PTS: m x (n-m) sized matrix --  --PST: (n-m) x m sized matrix --
makePTSandPST :: [Edge] -> Digraph -> Matrix Double -> (Matrix Double, Matrix Double)
makePTSandPST uVerts vs mat = (pts, pst)
  where vertices = nub $ concat $ map fstsnd vs  -- decomposed ant path
        idx      = accumIndex $ map nm $ map (\z -> elemIndex z vertices) uVerts -- for each uVert, get index.
        m        = (length uVerts) - 1
        nonidx   = (\\) [0..m] idx 
        prePST   = extractRows idx mat
        pst      = extractCols nonidx prePST
        prePTS   = extractRows nonidx mat
        pts      = extractCols idx prePTS

-- m x m sized matrix --
makePTT :: [Edge] -> Digraph -> Matrix Double -> Matrix Double
makePTT uVerts vs mat = ptt
  where vertices = nub $ concat $ map fstsnd vs  -- decomposed ant path
        idx      = filter (\w -> w > -1) $ map nm $ map (\z -> elemIndex z vertices) uVerts -- for each uVert, get index.
        m        = (length uVerts) - 1
        nonidx   = (\\) [1..m] idx 
        prePTT   = extractRows nonidx mat
        ptt      = extractCols nonidx prePTT

-- The absorbing model
-- G represents the probability of information arriving from an emitter.
scoreGraph :: State -> Ant -> ScoreData -> Double
scoreGraph s ant (uedges, mat) = x+y
  where (pts,pst) = makePTSandPST uedges (second ant) mat
        ptt = makePTT uedges (second ant) mat 
        i   = ident (cols ptt)
        z   = (inv (i - ptt))
        g   =  z<> pts
        h   = pst <> z
        x   = fromIntegral $ length $ filter (\x -> x > (tx s)) $ toList $ flatten h :: Double
        y   = fromIntegral $ length $ filter (\x -> x > (rx s)) $ toList $ flatten g :: Double        

-- The absorbing model
-- G represents the probability of information arriving from an emitter.
getScoreGraphMatrix :: Ant -> ScoreData -> Matrix Double
getScoreGraphMatrix ant (uedges, mat) = h
  where (pts,pst) = makePTSandPST uedges (second ant) mat
        ptt = makePTT uedges (second ant) mat 
        i   = ident (cols ptt)
        z   = (inv (i - ptt))
        g   =  z<> pts
        h   = pst <> z
        

score :: State -> Ant -> ScoreData -> Ant
score s ant scoredata
      | (hasDeadEdge ant) = ((node ant), (solution ant), 0.0)  -- trace ("dead edge")
      | otherwise = ((node ant), (solution ant), (scoreGraph s ant scoredata))

hasDeadEdge :: Ant -> Bool
hasDeadEdge (n1, es, sc) 
            | elem deadDedge es = True
            | otherwise = False 

-------------------------------------------------------------------------------------------------------------------------
-- END OF LINE --
-------------------------------------------------------------------------------------------------------------------------

-- ************************************************************************************************************************
