
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

module Biobase.Turner.Import.Turner where

import           Control.Arrow
import           Data.ByteString.Char8 (ByteString)
import           Data.ByteString.Lex.Fractional
import           Data.Char
import           Data.Default
import           Data.List (groupBy)
import           Data.List.Split
import           Data.Map.Strict (Map)
import           Data.Maybe (fromJust)
import qualified Data.ByteString.Char8 as BS
import qualified Data.List as L
import qualified Data.Map.Strict as M
import qualified Data.Vector.Unboxed as VU
import           System.FilePath.Posix

import           Biobase.Primary
import           Biobase.Primary.Letter
import           Biobase.Primary.Nuc.RNA
import           Biobase.Secondary
import           Biobase.Secondary.Convert
import           Biobase.Types.Energy
import           Data.PrimitiveArray hiding (C, map)
import qualified Biobase.Primary.Nuc.DNA as DNA

import           Biobase.Turner.Types



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
  multiMM'    <- blockFromFile $ fp </> prefix ++ "tstackm" <.> suffix
  imisc'      <- miscFromFile  $ fp </> prefix ++ "miscloop" <.> suffix
  extMM'      <- blockFromFile $ fp </> prefix ++ "tstack" <.> suffix
  coaxial'    <- blockFromFile $ fp </> prefix ++ "coaxial" <.> suffix
  cstack'     <- blockFromFile $ fp </> prefix ++ "coaxstack" <.> suffix
  tstack'     <- blockFromFile $ fp </> prefix ++ "tstackcoax" <.> suffix
  return Turner2004Model
    { _stack              = fromAssocs minBPP  maxBPP   def $ L.zip keysPP  stack'
    , _dangle3            = fromAssocs minBPB  maxBPB   def $ L.zip keysPB  dangle3'
    , _dangle5            = fromAssocs minBPB  maxBPB   def $ L.zip keysPB  dangle5'
    , _hairpinL           = VU.fromList $ def : hairpinL' -- fromAssocs (Z:.0) (Z:.30) def $ L.zip d1_30 hairpinL'
    , _hairpinMM          = fromAssocs minBPBB maxBPBB def $ L.zip keysPBB hairpinMM'
    , _hairpinLookup      = M.fromList $ hairpinLk3 ++ hairpinLk4 ++ hairpinLk6
    , _hairpinGGG         = DG . L.head $ imisc' !! 8
    , _hairpinCslope      = DG . L.head $ imisc' !! 9
    , _hairpinCintercept  = DG . L.head $ imisc' !! 10
    , _hairpinC3          = DG . L.head $ imisc' !! 11
    , _bulgeL             = VU.fromList $ def : bulgeL' -- fromAssocs (Z:.0)      (Z:.30)     def $ L.zip d1_30 bulgeL'
    , _bulgeSingleC       = DG . L.head $ imisc' !! 13
    , _iloop1x1           = fromAssocs minBPPBB   maxBPPBB   def $ L.zip keysPPBB   iloop1x1'
    , _iloop2x1           = fromAssocs minBPPBBB  maxBPPBBB  def $ L.zip keysPPBBB  iloop2x1'
    , _iloop2x2           = fromAssocs minBPPBBBB maxBPPBBBB def $ L.zip keysPPBBBBrna iloop2x2' -- (if (prefix == "" || suffix == "dh") then keysPPBBBBrna else keysPPBBBBdna) iloop2x2'
    , _iloopMM            = fromAssocs minBPBB    maxBPBB    def $ L.zip keysPBB    iloopMM'
    , _iloop2x3MM         = fromAssocs minBPBB    maxBPBB    def $ L.zip keysPBB    iloop2x3MM'
    , _iloop1xnMM         = fromAssocs minBPBB    maxBPBB    def $ L.zip keysPBB    iloop1xnMM'
    , _iloopL             = VU.fromList $ def : iloopL' -- fromAssocs (Z:.0)    (Z:.30)   def $ L.zip d1_30      iloopL'
    , _multiMM            = fromAssocs minBPBB    maxBPBB    def $ L.zip keysPBB    multiMM'
    , _ninio              = DG . L.head $ imisc' !! 2
    , _maxNinio           = DG . L.head $ imisc' !! 1
    , _multiOffset        = DG $ (imisc' !! 3) !! 0
    , _multiNuc           = DG $ (imisc' !! 3) !! 1
    , _multiHelix         = DG $ (imisc' !! 3) !! 2
    , _multiAsym          = DG . L.head $ imisc' !! 5
    , _multiStrain        = DG . L.head $ imisc' !! 6
    , _exteriorMM         = fromAssocs minBPBB maxBPBB def $ L.zip keysPBB extMM'
    , _coaxial            = fromAssocs minBPP  maxBPP  def $ L.zip keysPP  coaxial'
    , _coaxStack          = fromAssocs minBPBB maxBPBB def $ L.zip keysPBB cstack'
    , _tStackCoax         = fromAssocs minBPBB maxBPBB def $ L.zip keysPBB tstack'
    , _largeLoop          = DG . L.head $ imisc' !! 0
    , _termAU             = DG . L.head $ imisc' !! 7
    , _intermolecularInit = DG . L.head $ imisc' !! 12
    }

d1_30 = L.map (Z:.) [1..30]

-- note that the unorder in the indices is intentional!

keysPP     = [ Z:.basepairConvert (k1,k2):.basepairConvert(k4,k3)
             | k1 <- acgu, k3 <- acgu, k2 <- acgu, k4 <- acgu]
keysPB     = [ Z:.basepairConvert(k1,k2):.k3
             | k1 <- acgu, k2 <- acgu, k3 <- acgu]
keysPBB    = [ Z:.basepairConvert (k1,k2):.k3:.k4
             | k1 <- acgu, k3 <- acgu, k2 <- acgu, k4 <- acgu]
keysPPBB   = [ Z:.basepairConvert(k1,k2):.basepairConvert(k4,k3):.k5:.k6
             | (k1,k2) <- plist11, k5 <- acgu, (k3,k4) <- plist11, k6 <- acgu]
keysPPBBB  = [ Z:.basepairConvert(k1,k2):.basepairConvert(k4,k3):.k5:.k6:.k7
             | (k1,k2) <- plist11, k6 <- acgu, k5 <- acgu, (k3,k4) <- plist11, k7 <- acgu]
keysPPBBBBrna = [ Z:.basepairConvert(k1,k2):.basepairConvert(k4,k3):.k5:.k6:.k7:.k8
                | (k1,k2) <- plist22rna, (k3,k4) <- plist22rna, k5 <- acgu, k8 <- acgu, k6 <- acgu, k7 <- acgu]
keysPPBBBBdna = [ Z:.basepairConvert(k1,k2):.basepairConvert(k4,k3):.k5:.k6:.k7:.k8
                | (k1,k2) <- plist22dna, (k3,k4) <- plist22dna, k5 <- acgu, k8 <- acgu, k6 <- acgu, k7 <- acgu]

fullps = [basepairConvert a b | a <- acgu, b <- acgu]

plist11 = [(A,U),(C,G),(G,C),(U,A),(G,U),(U,G)]
plist22rna = [(A,U),(C,G),(G,C),(G,U),(U,A),(U,G)]
plist22dna = [(DNA.A,DNA.T),(DNA.C,DNA.G),(DNA.G,DNA.C),(DNA.T,DNA.A),(DNA.G,DNA.T),(DNA.T,DNA.G)]

-- * Conduit stuff

-- | extract values. "." - values are extracted as > 100k

values :: ByteString -> [DeltaGibbs]
values xs
  | BS.null ys
    = []
  | "." `BS.isPrefixOf` ys
    = def : values (BS.drop 1 ys)
  | Just (d,zs) <- readSigned readDecimal ys
    = DG d : values zs
  where ys = BS.dropWhile isSpace xs

-- | Iteratee to parse tabulated loops (hairpins).

parseTabulated :: ByteString -> [(ByteString,DeltaGibbs)]
parseTabulated = map f . filter (not . BS.all isSpace) . drop 2 . BS.lines
  where f x
          | Just (d,_) <- readSigned readDecimal v = (k,DG d)
          | otherwise = error $ "tabulated: <" ++ BS.unpack x ++ ">"
          where (k,v) = second (BS.dropWhile isSpace) . BS.span (not . isSpace) . BS.dropWhile isSpace $ x
{-
parseTabulated :: Monad m => Sink ByteString m [(ByteString,DeltaGibbs)]
parseTabulated = C.lines =$ CL.filter (not . BS.all isSpace) =$ g where
  g = do
    CL.drop 2
    xs <- CL.map f =$ consume
    return xs
  f x
    | Just (d,_) <- readDouble v = (k,DG d)
    | otherwise = error $ "tabulated: <" ++ BS.unpack x ++ ">"
    where (k,v) = second (BS.dropWhile isSpace) . BS.span (not . isSpace) . BS.dropWhile isSpace $ x
-}

-- | Convenience function

blockFromFile :: FilePath -> IO [DeltaGibbs]
blockFromFile = fmap parseBlocks . BS.readFile
{-
blockFromFile fp = do
  xs <- runResourceT $ sourceFile fp $$ parseBlocks =$ consume
  if (allEq $ L.map L.length xs)
    then return $ L.concat xs
    else error $ "in file: " ++ fp ++ " we have unequal line lengths"
-}

-- | Transform input stream into list of list of doubles

parseBlocks :: ByteString -> [DeltaGibbs]
parseBlocks = concat . filter (not . null) . map f . BS.lines
  where
{-
parseBlocks :: Monad m => Conduit ByteString m [DeltaGibbs]
parseBlocks = C.lines =$= CL.map f =$= CL.filter (not . L.null) where
-}
  f :: ByteString -> [DeltaGibbs]
  f x
    | "5'" `BS.isPrefixOf` y = []
    | "3'" `BS.isPrefixOf` y = []
    | "." `BS.isPrefixOf`  y = values y
    | Just (d,xs) <- readSigned readDecimal y = values y
    | otherwise = [] -- error $ BS.unpack x
    where y = BS.dropWhile isSpace x


-- | Parses the miscloop table
--
-- NOTE extra brownie points for miscloop.dat for providing data in a form that
-- does not conform to normal number encoding.

parseMiscLoop :: ByteString -> [[Double]]
parseMiscLoop = map f . drop 1 . groupBy (\x y -> not $ BS.null y) . BS.lines
  where f = map readD . BS.words . last
{-
parseMiscLoop :: Monad m => Sink ByteString m [[Double]]
parseMiscLoop = C.lines =$ CL.groupBy (\x y -> not $ BS.null y) =$ f where
  f = do
    CL.drop 1
    xs <- consume
    return . L.map (L.map readD . BS.words . L.last) $ xs
-}

-- | Parses stupidly encoded values like ".6" and "-.0".

readD :: ByteString -> Double
readD xs
  | BS.null xs              = error "readD: null"
  | BS.head xs == '.'       = readD $ BS.cons '0' xs
  | "-." `BS.isPrefixOf` xs = readD $ "-0." `BS.append` BS.drop 2 xs
  | Just (d,_) <- readSigned readDecimal xs = d
  | otherwise = error $ BS.unpack xs

-- |

miscFromFile :: FilePath -> IO [[Double]]
--miscFromFile fp = runResourceT $ sourceFile fp $$ parseMiscLoop
miscFromFile = fmap parseMiscLoop . BS.readFile

-- |

tabFromFile :: FilePath -> IO [(Primary RNA,DeltaGibbs)]
-- tabFromFile fp = fmap (L.map (first mkPrimary)) . runResourceT $ sourceFile fp $$ parseTabulated
tabFromFile = fmap (L.map (first primary)) . fmap parseTabulated . BS.readFile

allEq [] = True
allEq (x:xs) = L.all (==x) xs

type Prefix = FilePath
type Suffix = FilePath

