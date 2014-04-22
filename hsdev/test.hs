import Sparse
import Data.Maybe
import Data.IntMap.Strict as IM

subsetSM :: [Int] -> SM a -> [SV a]
subsetSM ks (SM (_, m)) = catMaybes $ Prelude.map (\x -> (x `IM.lookup` m)) ks

subsetSV :: [Int] -> SV a -> SV a
subsetSV [] (SV v) = (SV v)
subsetSV (x:ks) (SV v) = subsetSV ks (SV (IM.delete x v))


removeChosenK :: Int -> Int -> [Int] -> [Double] -> ([Int],[Double]) -> ([Int],[Double])
removeChosenK i j (ki:ks) (pi:ps) (x,y)  
  | j == 0 = (ks, ps)
  | i == j = (x ++ ks, y ++ ps)
  | otherwise = removeChosenK (i+1) j ks ps (x++[ki], y++[pi])


diffList :: [Int] -> [Int] -> [Int] -> [Int]
diffList [] ys zs = zs ++ ys
diffList _ [] zs = zs
diffList ks (y:ys) zs
  | elem y ks = diffList ks ys zs
  | otherwise = diffList ks ys (y:zs)

