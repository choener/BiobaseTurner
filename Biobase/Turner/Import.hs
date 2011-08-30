{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedStrings #-}

-- | Turner file parser. Returns a Turner2004 data structure. Requires an
-- annoying amount of boilerplate.
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
import Data.Maybe (fromJust)

import Biobase.Primary
import Biobase.Secondary
import Data.PrimitiveArray

import Biobase.Turner



-- *

-- | Given a directory, fill in the 'Turner2004' data structure

fromDir :: FilePath -> Prefix -> Suffix -> IO Turner2004
fromDir fp prefix suffix = do
  stack'     <- blockFromFile $ fp </> prefix ++ "stack"    <.> suffix
  dangle'    <- blockFromFile $ fp </> prefix ++ "dangle"   <.> suffix
  loop'      <- blockFromFile $ fp </> prefix ++ "loop"     <.> suffix
  hairpinMM' <- blockFromFile $ fp </> prefix ++ "tstackh"  <.> suffix
  hairpinLk3 <- tabFromFile   $ fp </> prefix ++ "triloop"  <.> suffix
  hairpinLk4 <- tabFromFile   $ fp </> prefix ++ "tloop"    <.> suffix
  hairpinLk6 <- tabFromFile   $ fp </> prefix ++ "hexaloop" <.> suffix
  let (dangle3',dangle5') = L.splitAt (L.length dangle' `div` 2) dangle'
  let (_:iloopL':bulgeL':hairpinL':[]) = L.transpose $ splitEvery 4 loop'
  iloop1x1'   <- blockFromFile $ fp </> prefix ++ "int11" <.> suffix
  iloop2x1'   <- blockFromFile $ fp </> prefix ++ "int21" <.> suffix
  iloop2x2'   <- blockFromFile $ fp </> prefix ++ "int22" <.> suffix
  iloopMM'    <- blockFromFile $ fp </> prefix ++ "tstacki" <.> suffix
  iloop2x3MM' <- blockFromFile $ fp </> prefix ++ "tstacki23" <.> suffix
  iloop1xnMM' <- blockFromFile $ fp </> prefix ++ "tstacki1n" <.> suffix
  multiMM'    <- blockFromFile $ fp </> prefix ++ "tstacki" <.> suffix
  imisc'      <- miscFromFile  $ fp </> prefix ++ "miscloop" <.> suffix
  extMM'      <- blockFromFile $ fp </> prefix ++ "tstack" <.> suffix
  coaxial'    <- blockFromFile $ fp </> prefix ++ "coaxial" <.> suffix
  cstack'     <- blockFromFile $ fp </> prefix ++ "coaxstack" <.> suffix
  tstack'     <- blockFromFile $ fp </> prefix ++ "tstackcoax" <.> suffix
  return Turner2004
    { stack         = fromAssocs minPP  maxPP  infE $ L.zip keysPP  stack'
    , dangle3       = fromAssocs minPB  maxPB  infE $ L.zip keysPB  dangle3'
    , dangle5       = fromAssocs minPB  maxPB  infE $ L.zip keysPB  dangle5'
    , hairpinL      = fromAssocs 0      30     infE $ L.zip [1..30] hairpinL'
    , hairpinMM     = fromAssocs minPBB maxPBB infE $ L.zip keysPBB hairpinMM'
    , hairpinLookup = M.fromList $ hairpinLk3 ++ hairpinLk4 ++ hairpinLk6
    , hairpinGGG    = L.head $ imisc' !! 8
    , hairpinCslope = L.head $ imisc' !! 9
    , hairpinCintercept = L.head $ imisc' !! 10
    , hairpinC3     = L.head $ imisc' !! 11
    , bulgeL        = fromAssocs 0      30     infE $ L.zip [1..30] bulgeL'
    , bulgeSingleC  = L.head $ imisc' !! 13
    , iloop1x1      = fromAssocs minPPBB   maxPPBB   infE $ L.zip keysPPBB   iloop1x1'
    , iloop2x1      = fromAssocs minPPBBB  maxPPBBB  infE $ L.zip keysPPBBB  iloop2x1'
    , iloop2x2      = fromAssocs minPPBBBB maxPPBBBB infE $ L.zip (if (prefix == "" || suffix == "dh") then keysPPBBBBrna else keysPPBBBBdna) iloop2x2'
    , iloopMM       = fromAssocs minPBB    maxPBB    infE $ L.zip keysPBB    iloopMM'
    , iloop2x3MM    = fromAssocs minPBB    maxPBB    infE $ L.zip keysPBB    iloop2x3MM'
    , iloop1xnMM    = fromAssocs minPBB    maxPBB    infE $ L.zip keysPBB    iloop1xnMM'
    , iloopL        = fromAssocs 0      30     infE $ L.zip [1..30] iloopL'
    , multiMM       = fromAssocs minPBB    maxPBB    infE $ L.zip keysPBB    multiMM'
    , ninio = L.head $ imisc' !! 2
    , maxNinio = L.head $ imisc' !! 1
    , multiOffset = (imisc' !! 3) !! 0
    , multiNuc = (imisc' !! 3) !! 1
    , multiHelix = (imisc' !! 3) !! 2
    , multiAsym = L.head $ imisc' !! 5
    , multiStrain = L.head $ imisc' !! 6
    , extMM      = fromAssocs minPBB maxPBB infE $ L.zip keysPBB extMM'
    , coaxial    = fromAssocs minPP  maxPP  infE $ L.zip keysPP  coaxial'
    , coaxStack  = fromAssocs minPBB maxPBB infE $ L.zip keysPBB cstack'
    , tStackCoax = fromAssocs minPBB maxPBB infE $ L.zip keysPBB tstack'
    , largeLoop = L.head $ imisc' !! 0
    , termAU = L.head $ imisc' !! 7
    , intermolecularInit = L.head $ imisc' !! 12
    }

minPP     = (minP,minP)
maxPP     = (maxP,maxP)
minP      = (nN,nN)
maxP      = (nU,nU)
minPB     = (minP,nN)
maxPB     = (maxP,nU)
minPBB    = (minP,nN,nN)
maxPBB    = (maxP,nU,nU)
minPPBB   = (minP,minP,(nN,nN))
maxPPBB   = (maxP,maxP,(nU,nU))
minPPBBB  = (minP,minP,(nN,nN,nN))
maxPPBBB  = (maxP,maxP,(nU,nU,nU))
minPPBBBB = (minP,minP,(nN,nN,nN,nN))
maxPPBBBB = (maxP,maxP,(nU,nU,nU,nU))

keysPP     = [((k1,k2),(k4,k3)) | k1 <- acgu, k3 <- acgu, k2 <- acgu, k4 <- acgu]
keysPB     = [((k1,k2),k3) | k1 <- acgu, k2 <- acgu, k3 <- acgu]
keysPBB    = [ ((k1,k2),k3,k4)
             | k1 <- acgu, k3 <- acgu, k2 <- acgu, k4 <- acgu]
keysPPBB   = [ ((k1,k2),(k4,k3),(k5,k6))
             | (k1,k2) <- plist11, k5 <- acgu, (k3,k4) <- plist11, k6 <- acgu]
keysPPBBB  = [ ((k1,k2),(k4,k3),(k5,k6,k7))
             | (k1,k2) <- plist11, k6 <- acgu, k5 <- acgu, (k3,k4) <- plist11, k7 <- acgu]
keysPPBBBBrna = [ ((k1,k2),(k4,k3),(k5,k6,k7,k8))
                | (k1,k2) <- plist22rna, (k3,k4) <- plist22rna, k5 <- acgu, k8 <- acgu, k6 <- acgu, k7 <- acgu]
keysPPBBBBdna = [ ((k1,k2),(k4,k3),(k5,k6,k7,k8))
                | (k1,k2) <- plist22dna, (k3,k4) <- plist22dna, k5 <- acgu, k8 <- acgu, k6 <- acgu, k7 <- acgu]

plist11 = [(nA,nU),(nC,nG),(nG,nC),(nU,nA),(nG,nU),(nU,nG)]
plist22rna = [(nA,nU),(nC,nG),(nG,nC),(nG,nU),(nU,nA),(nU,nG)]
plist22dna = [(nA,nT),(nC,nG),(nG,nC),(nT,nA),(nG,nT),(nT,nG)]

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

iTabulated :: (Functor m, Monad m) => Iteratee ByteString m [(ByteString,Double)]
iTabulated = joinI $ enumLinesBS ><> I.filter (BS.all isSpace) $ g where
  g = do
    I.drop 2
    joinI $ mapStream f stream2stream
  f x
    | Just (d,_) <- readDouble v = (k,d)
    | otherwise = error $ "tabulated: <" ++ BS.unpack x ++ ">"
    where (k,v) = second (BS.dropWhile isSpace) . BS.span (not . isSpace) . BS.dropWhile isSpace $ x

-- | Convenience function

blockFromFile :: FilePath -> IO [Double]
blockFromFile fp = do
  i <- enumFile 8192 fp . joinI $ eneeBlocks stream2list
  xs <- run i
  if (allEq $ L.map L.length xs)
    then return $ L.concat xs
    else error $ "in file: " ++ fp ++ " we have unequal line lengths"

-- | Parses the miscloop table
--
-- NOTO extra brownie points for miscloop.dat for providing data in a form that
-- does not conform to normal number encoding.

iMiscLoop :: (Functor m, Monad m) => Iteratee ByteString m [[Double]]
iMiscLoop = joinI $ enumLinesBS ><> I.groupBy (\x y -> not $ BS.null y) $ f where
  f = do
    I.drop 1
    xs <- fmap (L.map (L.map (readD) . BS.words . L.last)) $ stream2list
    return xs

-- | Parses stupidly encoded values like ".6" and "-.0".

readD :: ByteString -> Double
readD xs
  | BS.null xs              = error "readD: null"
  | BS.head xs == '.'       = readD $ BS.cons '0' xs
  | "-." `BS.isPrefixOf` xs = readD $ "-0." `BS.append` BS.drop 2 xs
  | Just (d,_) <- readDouble xs = d
  | otherwise = error $ BS.unpack xs

-- |

miscFromFile :: FilePath -> IO [[Double]]
miscFromFile fp = run =<< enumFile 8192 fp iMiscLoop

-- | 

tabFromFile :: FilePath -> IO [(ByteString,Double)]
tabFromFile fp = run =<< enumFile 8192 fp iTabulated

allEq [] = True
allEq (x:xs) = L.all (==x) xs

{-
test1 = blockFromFile "/home/choener/Documents/Workdata/TurnerRNA2004/RNA/loop.dat"
test2 = fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
-}

type Prefix = FilePath
type Suffix = FilePath

