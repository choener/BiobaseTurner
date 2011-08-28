{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedStrings #-}

-- | Turner file parser. Returns a Turner2004 data structure. We store data in
-- the same way it is stored in the ViennaRNA package. Pairs are tuples
-- however.
--
-- How is 'stack' data stored:
--
-- AX
-- UY   ->  ((A,U),(Y,X))
--
-- How 'iloop1x1' is stored:
--
--  X
-- A G
-- U C  -> ((A,U),(C,G),X,Y)
--  Y
--
--  Now 'iloop1x2' is stored:
--
--  X
-- A  G
-- U  C  -> ((A,U),(C,G),X,C,Y), single (X) first, then 5' to 3'
--  YC
--
--  'iloop2x2' is stored:
--
--  XY
-- A  G
-- U  C  -> ((A,U),(C,G),X,Y,y,x), X-->Y then x<--y
--  xy
--
--TODO not sure if dangle3/dangle5 are correctly split or if they should switch

module Biobase.Turner.Import where

import Control.Arrow
import Data.ByteString.Char8 as BS
import Data.ByteString.Lex.Double
import Data.Char
import Data.Iteratee as I
import Data.Iteratee.Char as I
import Data.Iteratee.IO as I
import Data.List.Split
import Data.Map as M
import qualified Data.List as L
import System.FilePath.Posix

import Biobase.Primary
import Biobase.Secondary
import Data.PrimitiveArray

import Biobase.Turner



-- *

-- | Given a directory, fill in the 'Turner2004' data structure

fromDir :: FilePath -> FilePath -> FilePath -> IO Turner2004
fromDir fp prefix suffix = do
  stack'     <- blockFromFile $ fp </> prefix ++ "stack"    <.> suffix
  dangle'    <- blockFromFile $ fp </> prefix ++ "dangle"   <.> suffix
  loop'      <- blockFromFile $ fp </> prefix ++ "loop"     <.> suffix
  hairpinMM' <- blockFromFile $ fp </> prefix ++ "tstackh"  <.> suffix
  hairpinLk3 <- tabFromFile   $ fp </> prefix ++ "triloop"  <.> suffix
  hairpinLk4 <- tabFromFile   $ fp </> prefix ++ "tloop"    <.> suffix
  hairpinLk6 <- tabFromFile   $ fp </> prefix ++ "hexaloop" <.> suffix
  let (dangle3',dangle5') = L.splitAt (L.length dangle' `div` 2) dangle'
  let (_:viloopl:vbulgel:hairpinL':[]) = L.transpose $ splitEvery 4 loop'
  return Turner2004
    { stack         = fromAssocs minPP  maxPP  infE $ L.zip keysPP  stack'
    , dangle3       = fromAssocs minPB  maxPB  infE $ L.zip keysPB  dangle3'
    , dangle5       = fromAssocs minPB  maxPB  infE $ L.zip keysPB  dangle5'
    , hairpinL      = fromAssocs 0      30     infE $ L.zip [1..30] hairpinL'
    , hairpinMM     = fromAssocs minPBB maxPBB infE $ L.zip keysPBB hairpinMM'
    , hairpinLookup = M.fromList $ hairpinLk3 ++ hairpinLk4 ++ hairpinLk6
    }

minPP  = (minP,minP)
maxPP  = (maxP,maxP)
minP   = (nN,nN)
maxP   = (nU,nU)
minPB  = (minP,nN)
maxPB  = (maxP,nU)
minPBB = (minP,nN,nN)
maxPBB = (maxP,nU,nU)

keysPP = [((k1,k2),(k4,k3)) | k1 <- acgu, k3 <- acgu, k2 <- acgu, k4 <- acgu]
keysPB = [((k1,k2),k3) | k1 <- acgu, k2 <- acgu, k3 <- acgu]
keysPBB = [((k1,k2),k3,k4) | k1 <- acgu, k3 <- acgu, k2 <- acgu, k4 <- acgu]

infE = 999999 :: Double

-- * Iteratee stuff

-- | Transform input stream into list of list of doubles

eneeBlocks :: (Functor m, Monad m) => Enumeratee ByteString [[Double]] m a
eneeBlocks = enumLinesBS ><> mapStream f ><> I.filter (not . L.null) where
  f x
    | "5'" `isPrefixOf` y = []
    | "3'" `isPrefixOf` y = []
    | "." `isPrefixOf`  y = values y
    | Just (d,xs) <- readDouble y = values y
    | otherwise = [] -- error $ BS.unpack x
    where y = BS.dropWhile isSpace x

-- | extract values. "." - values are extracted as > 100k

values :: ByteString -> [Double]
values xs
  | BS.null ys
    = []
  | "." `isPrefixOf` ys
    = infE : values (BS.drop 1 ys)
  | Just (d,zs) <- readDouble ys
    = d : values zs
  where ys = BS.dropWhile isSpace xs

-- | Iteratee to parse tabulated loops (hairpins).

iTabulated :: (Monad m) => Iteratee ByteString m [(ByteString,Double)]
iTabulated = joinI $ enumLinesBS g where
  g = do
    I.drop 2
    joinI $ mapStream f stream2stream
  f x
    | Just (d,_) <- readDouble v = (k,d)
    | otherwise = error $ BS.unpack x
    where (k,v) = second (BS.dropWhile isSpace) . BS.span (not . isSpace) . BS.dropWhile isSpace $ x

-- | Convenience function

blockFromFile :: FilePath -> IO [Double]
blockFromFile fp = do
  i <- enumFile 8192 fp . joinI $ eneeBlocks stream2list
  xs <- run i
  if (allEq $ L.map L.length xs)
    then return $ L.concat xs
    else error $ "in file: " ++ fp ++ " we have unequal line lengths"



tabFromFile :: FilePath -> IO [(ByteString,Double)]
tabFromFile fp = run =<< enumFile 8192 fp iTabulated

allEq [] = True
allEq (x:xs) = L.all (==x) xs

test1 = blockFromFile "/home/choener/Documents/Workdata/TurnerRNA2004/RNA/loop.dat"

