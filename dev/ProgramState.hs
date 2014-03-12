-------------------------------------------------------------------------------------------------------------------------
-- Influence Maximization Problem using hypercube-minmax ant optimization ----------------------------------------------
-- for problems on biological pathway graphs 
-------------------------------------------------------------------------------------------------------------------------

-- main = getArgs >>=  printResult . antopt . buildGraphs . buildState
--                                                             ***
-- Input:
-- String list from getArgs --

-- Output:
-- Program state record
-------------------------------------------------------------------------------------------------------------------------

module ProgramState 
(State
,initState
,showState
,parseValue
,readState
,fillState
,emptyState
,isEmptyState
,damp
,runs,ants,ct,tx,rx,local,evap,alph,beta
,bestEver,bestRest,bestIter, ptoKratio
,g,config,graphfile,randomSeed,c,k,runsDone
) where

import System.IO.Unsafe
import Data.Char (isSpace)
import Utilities
import System.Random
import Numeric.LinearAlgebra
import Data.List (nub)

data State = 
  State {    runs :: Int,  -- the number of restarts
             ants :: Int,  -- the number of ants
             ct :: Double, -- converge threshold
             tx :: Double, -- the transmit threshold
             rx :: Double, -- the receive threshold
             local :: Int, -- turn local optimization on?
             damp :: Double, -- the level of dampening on signal prop
             evap :: Double, -- the evaporation rate
             alph :: Double,  -- alpha in the ant op. algorithm
             beta :: Double,  -- beta in the ant op. algorithm
             bestEver :: [Int], -- best over all
             bestRest :: [Int], -- best since the last restart
             bestIter :: [Int], -- best this ant-iteration.
             g  :: StdGen, -- the random generator
             config :: String, -- the configuration file
             graphfile :: String, -- the graph data file
             randomSeed :: Int, -- the random seed -- 0 indicates a new generator.
             c :: Double, -- the convergence value
             k :: Double,  -- the number of nodes to select!
             runsDone :: Int,  -- the current ant
             ptoKratio :: Double -- the phenomone sum / k
        }

parseValue :: String -> String
parseValue s = parseValue' 0 s []

parseValue' :: Int -> String -> String -> String
parseValue' 1 [] ys = reverse ys
parseValue' flag (x:xs) ys 
  | flag == 0 && x /= ':' = parseValue' 0 xs ys
  | flag == 0 && x == ':' = parseValue' 1 xs ys
  | flag == 1 && (isSpace x) = parseValue' 1 xs ys 
  | flag == 1 && not (isSpace x) = parseValue' 1 xs (x:ys)

initState :: (StdGen, [String]) -> Maybe State
initState (g,[])       = Nothing
initState (g,(x:[]))   = Nothing
initState (g,(x:y:[])) = Just $ fillState g (readState x) x y
initState (g,(x:y:xs)) = Nothing

readState :: String -> [String] 
readState s = map parseValue $ lines $ unsafePerformIO $ readFile s  

fillState :: StdGen -> [String] -> String -> String -> State
fillState g xs y z
  | (xs !! 11) == "Random" = fillState1 g xs y z
  | otherwise = fillState2 xs y z

fillState1 :: StdGen -> [String] -> String -> String -> State
fillState1 gen xs y z = State {runs = read (xs !! 1) :: Int,
                          ants = read (xs !! 2) :: Int,
                          ct   = read (xs !! 3) :: Double, 
                          local = read (xs !! 4) :: Int,
                          evap = read (xs !! 5) :: Double,
                          damp = read (xs !! 6) :: Double,
                          alph  = read (xs !! 7) :: Double,
                          beta  = read (xs !! 8) :: Double,
                          tx = read (xs !! 9) :: Double,
                          rx = read (xs !! 10) :: Double,
                          bestEver = [] :: [Int],
                          bestRest = [] :: [Int],
                          bestIter = [] :: [Int],
                          g = gen,
                          randomSeed = 0,
                          config = y,
                          graphfile = z,
                          c = 1.0,
                          k = (read (xs !! 12) :: Double),
                          ptoKratio = 0.0, 
                          runsDone = 0
                     }

showState :: State -> Digraph -> String
showState s digraph = 
  (show (config s)) ++ "\t" ++ 
  (show (graphfile s)) ++ "\t" ++ 
  (show (randomSeed s)) ++ "\t" ++
  (show (k s)) ++ "\t" ++ 
  (show (runs s)) ++ "\t" ++ 
  (show (ants s)) ++ "\t" ++ 
  (show (ct s)) ++ "\t" ++
  (show (alph s)) ++ "\t" ++
  (show (beta s)) ++ "\t" ++
  (show (local s)) ++ "\t" ++ 
  (show (damp s)) ++ "\t" ++
  (show (evap s)) ++ "\t" ++ 
  (show (tx s)) ++ "\t" ++
  (show (rx s)) ++ "\t" ++ 
  (show (ptoKratio s)) ++ "\t" ++
  (show (bestEver s)) ++ "\n"

  
-- prettyPrintGraph :: Digraph -> String -> String
-- prettyPrintGraph [] str = str 
-- prettyPrintGraph (d:dig) str = prettyPrintGraph dig (((show d) ++ ":") ++ str)


fillState2 :: [String] -> String -> String -> State
fillState2 xs y z = State {runs = read (xs !! 1) :: Int,
                          ants = read (xs !! 2) :: Int,
                          ct   = read (xs !! 3) :: Double, 
                          local = read (xs !! 4) :: Int,
                          evap = read (xs !! 5) :: Double,
                          damp = read (xs !! 6) :: Double,
                          alph  = read (xs !! 7) :: Double,
                          beta  = read (xs !! 8) :: Double,
                          tx = read (xs !! 9) :: Double,
                          rx = read (xs !! 10) :: Double,
                          bestEver = [] :: [Int],
                          bestRest = [] :: [Int],
                          bestIter = [] :: [Int],
                          g = mkStdGen (read (xs !! 11) :: Int),
                          randomSeed = (read (xs !! 11) :: Int),
                          config = y,
                          graphfile = z,
                          c = 1.0,
                          k = (read (xs !! 12) :: Double),
                          runsDone = 0, 
                          ptoKratio = 0.0
                          }


isEmptyState :: State -> Bool
isEmptyState s 
  | (runs s == 0) && (ants s == 0) && (config s == "x") && (graphfile s == "y") = True
  | otherwise = False

emptyState :: State
emptyState = State {runs = 0 :: Int,
                    ants = 0 :: Int,
                    ct   = 0.0 :: Double, 
                    local = 0 :: Int,
                    evap = 0.0 :: Double,
                    damp = 0.0 :: Double,
                    alph  = 0.0 :: Double,
                    beta  = 0.0 :: Double,
                    tx = 0.0 :: Double,
                    rx = 0.0 :: Double,
                    bestEver = [] :: [Int],
                    bestRest = [] :: [Int] ,
                    bestIter = [] :: [Int],
                    g = mkStdGen 100,
                    config = "x",
                    graphfile = "y",
                    c = 0.0,
                    k = 0.0,
                    randomSeed=100,
                    runsDone = 0,
                    ptoKratio = 0.0
                   }
             


-------------------------------------------------------------------------------------------------------------------------
-- END OF LINE --
-------------------------------------------------------------------------------------------------------------------------

-- ************************************************************************************************************************

