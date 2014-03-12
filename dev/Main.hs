-------------------------------------------------------------------------------------------------------------------------
-- Influence Maximization Problem using hypercube-minmax ant optimization ----------------------------------------------
-- for problems on biological pathway graphs 
-------------------------------------------------------------------------------------------------------------------------

-- Input:
-- config file and directed pathways with edge weights --

-- Output:
-- The program state
-- approximate best subset of size K, with the greatest influence on graph, by diffusion models.
-------------------------------------------------------------------------------------------------------------------------


module Main
       (main,
        printResult) where
  
import System.Environment
import System.Random
import Utilities
import ProgramState
import Graph
--import AntOpt
--import InfoFlow
import Graph


printResult :: String -> IO()
printResult result = putStrLn result

main :: IO ()
main = newStdGen >>= \g -> getArgs >>= \args -> 
  (printResult . initGraphDataT . initState) (g,args)

--main :: IO ()
--main = newStdGen >>= \g -> getArgs >>= \args -> 
-- (printResult . optimize . initGraphData . initState) (g,args)
