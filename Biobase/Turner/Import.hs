{-# LANGUAGE PatternGuards #-}
{-# LANGUAGE TupleSections #-}
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
import Data.Array.Repa.Index
import Data.ByteString.Char8 as BS
import Data.ByteString.Lex.Double
import Data.Char
import Data.Conduit as C
import Data.Conduit.Binary as C
import Data.Conduit.List as CL
import Data.List.Split
import Data.Map as M
import Data.Maybe (fromJust)
import qualified Data.List as L
import System.FilePath.Posix

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
  let (_:iloopL':bulgeL':hairpinL':[]) = L.transpose $ chunksOf 4 loop'
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
  return Turner2004Model
    { _stack              = fromAssocs minPP  maxPP   infE $ L.zip keysPP  stack'
    , _dangle3            = fromAssocs minPB  maxPB   infE $ L.zip keysPB  dangle3'
    , _dangle5            = fromAssocs minPB  maxPB   infE $ L.zip keysPB  dangle5'
    , _hairpinL           = fromAssocs (Z:.0) (Z:.30) infE $ L.zip d1_30 hairpinL'
    , _hairpinMM          = fromAssocs minPBB maxPBB infE $ L.zip keysPBB hairpinMM'
    , _hairpinLookup      = M.fromList $ hairpinLk3 ++ hairpinLk4 ++ hairpinLk6
    , _hairpinGGG         = Energy . L.head $ imisc' !! 8
    , _hairpinCslope      = Energy . L.head $ imisc' !! 9
    , _hairpinCintercept  = Energy . L.head $ imisc' !! 10
    , _hairpinC3          = Energy . L.head $ imisc' !! 11
    , _bulgeL             = fromAssocs (Z:.0)      (Z:.30)     infE $ L.zip d1_30 bulgeL'
    , _bulgeSingleC       = Energy . L.head $ imisc' !! 13
    , _iloop1x1           = fromAssocs minPPBB   maxPPBB   infE $ L.zip keysPPBB   iloop1x1'
    , _iloop2x1           = fromAssocs minPPBBB  maxPPBBB  infE $ L.zip keysPPBBB  iloop2x1'
    , _iloop2x2           = fromAssocs minPPBBBB maxPPBBBB infE $ L.zip (if (prefix == "" || suffix == "dh") then keysPPBBBBrna else keysPPBBBBdna) iloop2x2'
    , _iloopMM            = fromAssocs minPBB    maxPBB    infE $ L.zip keysPBB    iloopMM'
    , _iloop2x3MM         = fromAssocs minPBB    maxPBB    infE $ L.zip keysPBB    iloop2x3MM'
    , _iloop1xnMM         = fromAssocs minPBB    maxPBB    infE $ L.zip keysPBB    iloop1xnMM'
    , _iloopL             = fromAssocs (Z:.0)    (Z:.30)   infE $ L.zip d1_30      iloopL'
    , _multiMM            = fromAssocs minPBB    maxPBB    infE $ L.zip keysPBB    multiMM'
    , _ninio              = Energy . L.head $ imisc' !! 2
    , _maxNinio           = Energy . L.head $ imisc' !! 1
    , _multiOffset        = Energy $ (imisc' !! 3) !! 0
    , _multiNuc           = Energy $ (imisc' !! 3) !! 1
    , _multiHelix         = Energy $ (imisc' !! 3) !! 2
    , _multiAsym          = Energy . L.head $ imisc' !! 5
    , _multiStrain        = Energy . L.head $ imisc' !! 6
    , _extMM              = fromAssocs minPBB maxPBB infE $ L.zip keysPBB extMM'
    , _coaxial            = fromAssocs minPP  maxPP  infE $ L.zip keysPP  coaxial'
    , _coaxStack          = fromAssocs minPBB maxPBB infE $ L.zip keysPBB cstack'
    , _tStackCoax         = fromAssocs minPBB maxPBB infE $ L.zip keysPBB tstack'
    , _largeLoop          = Energy . L.head $ imisc' !! 0
    , _termAU             = Energy . L.head $ imisc' !! 7
    , _intermolecularInit = Energy . L.head $ imisc' !! 12
    }

minPP     = Z:.nN:.nN:.nN:.nN -- (minP,minP)
maxPP     = Z:.nU:.nU:.nU:.nU -- (maxP,maxP)
minP      = Z:.nN:.nN -- (nN,nN)
maxP      = Z:.nU:.nU -- (nU,nU)
minPB     = minP:.nN -- (minP,nN)
maxPB     = maxP:.nU -- (maxP,nU)
minPBB    = minPB:.nN -- (minP,nN,nN)
maxPBB    = maxPB:.nU -- (maxP,nU,nU)
minPPBB   = minPP:.nN:.nN -- (minP,minP,(nN,nN))
maxPPBB   = maxPP:.nU:.nU -- (maxP,maxP,(nU,nU))
minPPBBB  = minPPBB:.nN -- (minP,minP,(nN,nN,nN))
maxPPBBB  = maxPPBB:.nU -- (maxP,maxP,(nU,nU,nU))
minPPBBBB = minPPBBB:.nN -- (minP,minP,(nN,nN,nN,nN))
maxPPBBBB = maxPPBBB:.nU -- (maxP,maxP,(nU,nU,nU,nU))

d1_30 = L.map (Z:.) [1..30]

keysPP     = [{- ((k1,k2),(k4,k3)) -} Z:.k1:.k2:.k4:.k3 | k1 <- acgu, k3 <- acgu, k2 <- acgu, k4 <- acgu]
keysPB     = [{- ((k1,k2),k3) -} Z:.k1:.k2:.k3 | k1 <- acgu, k2 <- acgu, k3 <- acgu]
keysPBB    = [ {- ((k1,k2),k3,k4) -} Z:.k1:.k2:.k3:.k4
             | k1 <- acgu, k3 <- acgu, k2 <- acgu, k4 <- acgu]
keysPPBB   = [ {- ((k1,k2),(k4,k3),(k5,k6)) -} Z:.k1:.k2:.k4:.k3:.k5:.k6
             | (k1,k2) <- plist11, k5 <- acgu, (k3,k4) <- plist11, k6 <- acgu]
keysPPBBB  = [ {- ((k1,k2),(k4,k3),(k5,k6,k7)) -} Z:.k1:.k2:.k4:.k3:.k5:.k6:.k7
             | (k1,k2) <- plist11, k6 <- acgu, k5 <- acgu, (k3,k4) <- plist11, k7 <- acgu]
keysPPBBBBrna = [ {- ((k1,k2),(k4,k3),(k5,k6,k7,k8)) -} Z:.k1:.k2:.k4:.k3:.k5:.k6:.k7:.k8
                | (k1,k2) <- plist22rna, (k3,k4) <- plist22rna, k5 <- acgu, k8 <- acgu, k6 <- acgu, k7 <- acgu]
keysPPBBBBdna = [ {- ((k1,k2),(k4,k3),(k5,k6,k7,k8)) -} Z:.k1:.k2:.k4:.k3:.k5:.k6:.k7:.k8
                | (k1,k2) <- plist22dna, (k3,k4) <- plist22dna, k5 <- acgu, k8 <- acgu, k6 <- acgu, k7 <- acgu]

plist11 = [(nA,nU),(nC,nG),(nG,nC),(nU,nA),(nG,nU),(nU,nG)]
plist22rna = [(nA,nU),(nC,nG),(nG,nC),(nG,nU),(nU,nA),(nU,nG)]
plist22dna = [(nA,nT),(nC,nG),(nG,nC),(nT,nA),(nG,nT),(nT,nG)]

infE = Energy 999999

-- * Conduit stuff

-- | extract values. "." - values are extracted as > 100k

values :: ByteString -> [Energy]
values xs
  | BS.null ys
    = []
  | "." `isPrefixOf` ys
    = infE : values (BS.drop 1 ys)
  | Just (d,zs) <- readDouble ys
    = Energy d : values zs
  where ys = BS.dropWhile isSpace xs

-- | Iteratee to parse tabulated loops (hairpins).

parseTabulated :: Monad m => Sink ByteString m [(ByteString,Energy)]
parseTabulated = C.lines =$ CL.filter (not . BS.all isSpace) =$ g where
  g = do
    CL.drop 2
    xs <- CL.map f =$ consume
    return xs
  f x
    | Just (d,_) <- readDouble v = (k,Energy d)
    | otherwise = error $ "tabulated: <" ++ BS.unpack x ++ ">"
    where (k,v) = second (BS.dropWhile isSpace) . BS.span (not . isSpace) . BS.dropWhile isSpace $ x

-- | Convenience function

blockFromFile :: FilePath -> IO [Energy]
blockFromFile fp = do
  xs <- runResourceT $ sourceFile fp $$ parseBlocks =$ consume
  if (allEq $ L.map L.length xs)
    then return $ L.concat xs
    else error $ "in file: " ++ fp ++ " we have unequal line lengths"

-- | Transform input stream into list of list of doubles

parseBlocks :: Monad m => Conduit ByteString m [Energy]
parseBlocks = C.lines =$= CL.map f =$= CL.filter (not . L.null) where
  f x
    | "5'" `isPrefixOf` y = []
    | "3'" `isPrefixOf` y = []
    | "." `isPrefixOf`  y = values y
    | Just (d,xs) <- readDouble y = values y
    | otherwise = [] -- error $ BS.unpack x
    where y = BS.dropWhile isSpace x



-- | Parses the miscloop table
--
-- NOTE extra brownie points for miscloop.dat for providing data in a form that
-- does not conform to normal number encoding.

parseMiscLoop :: Monad m => Sink ByteString m [[Double]]
parseMiscLoop = C.lines =$ CL.groupBy (\x y -> not $ BS.null y) =$ f where
  f = do
    CL.drop 1
    xs <- consume
    return . L.map (L.map readD . BS.words . L.last) $ xs

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
miscFromFile fp = runResourceT $ sourceFile fp $$ parseMiscLoop

-- |

tabFromFile :: FilePath -> IO [(ByteString,Energy)]
tabFromFile fp = runResourceT $ sourceFile fp $$ parseTabulated

allEq [] = True
allEq (x:xs) = L.all (==x) xs

type Prefix = FilePath
type Suffix = FilePath

