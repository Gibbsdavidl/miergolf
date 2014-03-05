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
import AntOpt
import InfoFlow
import Graph

header :: String
header =   
  "Config\t" ++ 
  "Graph\t" ++ 
  "Nodes\t" ++
  "Edges\t" ++
  "Random_Seed\t" ++
  "K\t" ++ 
  "Runs\t" ++ 
  "Ants\t" ++ 
  "Convergence_Threshold\t" ++
  "Alpha\t" ++
  "Beta\t" ++
  "Local_Optim\t" ++ 
  "Dampening\t" ++
  "Evaporation\t" ++ 
  "TX_Threshold\t" ++
  "RX_Threshold\t" ++ 
  "PtoKRatio\t" ++
  "Score\t" ++
  "BestSolution" ++ "\n"


printResult :: String -> IO()
printResult result = putStrLn $ header ++ "\n" ++ result

main :: IO ()
main = newStdGen >>= \g -> getArgs >>= \args -> 
  (printResult . optimize . initGraphData . initState) (g,args)
