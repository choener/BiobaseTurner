
-- | Import a @ViennaRNA 2004@ energy parameter file.

module Biobase.Turner.Import.Vienna where

import Control.Applicative ( (<|>), Alternative )
import Control.Carrier.Lift
import Control.Carrier.Throw.Either
import Control.Effect.Throw
import Control.Monad (guard, MonadPlus)
import Control.Monad.IO.Class (MonadIO, liftIO)
import Data.ByteString (ByteString)
import Data.Char (isSpace)
import Data.List (foldl')
import Data.Map.Strict (Map)
import Data.PrimitiveArray as PA hiding (map,fromList)
import Data.Semigroup ((<>))
import Data.Text (Text,pack)
import Data.Vector.Unboxed (Vector,fromList)
import qualified Data.Map.Strict as M
import Text.Parser.Token.Style
import Text.Printf
import Text.Trifecta as TT
import Text.Trifecta.Delta as TD

import Algebra.Structure.Semiring
import Biobase.Primary.Letter
import Biobase.Primary.Nuc.RNA
import Biobase.Types.BioSequence
import Biobase.Types.Energy

import Biobase.Turner.Types



type R = Letter RNA ()
type Vienna2004 = Turner2004 Vector Vector R DDG

newtype P a = P { runP :: Parser a }
  deriving newtype ( Functor, Applicative, Monad, MonadPlus, Alternative, Parsing, CharParsing, DeltaParsing )

instance TokenParsing P where
  someSpace = buildSomeSpaceParser (skipSome (satisfy isSpace)) javaCommentStyle

-- | Returns two model structures, one with energies, one with enthalpies.
--
-- TODO Use @data ErrorString (s::String) | ErrorPrint (p :: IO())@ to allow for pretty-printing
-- errors

fromFile :: (MonadIO m, Has (Throw String) sig m) => FilePath -> m (Vienna2004, Vienna2004)
fromFile fp = TT.parseFromFileEx pVienna fp >>= \case
  Success r -> return r
  Failure e -> throwError ("Biobase.Turner.Import.Vienna.fromFile" ++ show e)

-- | Parses a given bytestring. Errors are rendered to a simple string for now. Mostly for use with
-- embedding.

fromByteString :: ByteString -> Either String (Vienna2004,Vienna2004)
fromByteString bs = case TT.parseByteString pVienna mempty bs of
  Success r -> Right r
  Failure e -> Left $ show e

-- | High-level parser that also detect correct conclusion of parsing.

pVienna :: Parser (Vienna2004,Vienna2004)
pVienna = runP $ v2Header *> spaces *> energiesAndEnthalpies <* v2End <* eof
  where v2Header = text "## RNAfold parameter file v2.0" <?> ".par Header"
        v2End = text "#END" <* newline
        energiesAndEnthalpies = do
          whiteSpace
          s  <- Stack <$> blockPP "# stack"
          sE <- Stack <$> blockPP "# stack_enthalpies"
          mmh  <- blockPNN "# mismatch_hairpin"
          mmhE <- blockPNN "# mismatch_hairpin_enthalpies"
          mmi  <- blockPNN "# mismatch_interior"
          mmiE <- blockPNN "# mismatch_interior_enthalpies"
          mmi1n  <- blockPNN "# mismatch_interior_1n"
          mmi1nE <- blockPNN "# mismatch_interior_1n_enthalpies"
          mmi23  <- blockPNN "# mismatch_interior_23"
          mmi23E <- blockPNN "# mismatch_interior_23_enthalpies"
          mmm  <- blockPNN "# mismatch_multi"
          mmmE <- blockPNN "# mismatch_multi_enthalpies"
          mme  <- blockPNN "# mismatch_exterior"
          mmeE <- blockPNN "# mismatch_exterior_enthalpies"
          d5  <- blockPN "# dangle5"
          d5E <- blockPN "# dangle5_enthalpies"
          d3  <- blockPN "# dangle3"
          d3E <- blockPN "# dangle3_enthalpies"
          i11  <- blockPPNN "# int11"
          i11E <- blockPPNN "# int11_enthalpies"
          i21  <- blockPPNNN "# int21"
          i21E <- blockPPNNN "# int21_enthalpies"
          i22  <- blockPPNNNN "# int22"
          i22E <- blockPPNNNN "# int22_enthalpies"
          hpl  <- blockL "# hairpin"
          hplE <- blockL "# hairpin_enthalpies"
          bul  <- blockL "# bulge"
          bulE <- blockL "# bulge_enthalpies"
          int  <- blockL "# interior"
          intE <- blockL "# interior_enthalpies"
          -- NOTE: combined entropy/enthalpy parameters
          (ml,mlE) <- string "# ML_params" *> someSpace *> ((,) <$> pmlLoop <*> pmlLoop)
          -- TODO: ignored ninio, misc
          __ninio <- string "# NINIO" *> someSpace *> some ddg
          __misc <- string "# Misc" *> someSpace *> some ddg
          (hexloop,hexloopE) <- loopMaps "# Hexaloops"
          (tetraloop,tetraloopE) <- loopMaps "# Tetraloops"
          (triloop,triloopE) <- loopMaps "# Triloops"
          -- NOTE end penalties are derived
          let
            penalty = __misc !! 2
            end = PA.fromAssocs maxBound one $ map (,penalty) [ (Z:.A:.U), (Z:.U:.A), (Z:.G:.U), (Z:.U:.G) ]
            endE = PA.fromAssocs maxBound one [] -- TODO figure out the enthalpies
          let
            hairpin = Hairpin
              { _hairpinLength = hpl
              , _hairpinMM = mmh
              , _hairpinLookup = hexloop <> tetraloop <> triloop
              , _hairpinPenalty = end -- NOTE same penalty as in exterior loops
--              , _hairpinGGG = ddgDbg
--              , _hairpinCslope = ddgDbg
--              , _hairpinCintercept = ddgDbg
--              , _hairpinC3 = ddgDbg
              , _largeLoop = ddgDbg
              }
            intloop = IntLoop
              { _intLoop1x1 = i11
              , _intLoop1x2 = i21
              , _intLoop2x2 = i22
              , _intLoop1xn = mmi1n
              , _intLoopMM  = mmi
              , _intLoopL = int
              , _bulgeL = bul
              }
            multi = Multi
              { _mismatchMulti = mmm
              }
            entropy = Turner2004
              { _hairpin = hairpin
              , _stack   = s
              , _intloop = intloop
              , _multi   = multi
              , _exterior = Exterior
                { _mismatchExterior = mme
                , _dangle5 = d5
                , _dangle3 = d3
                , _endPenalty = end
                }
              }
            hairpinE = Hairpin
              { _hairpinLength = hplE
              , _hairpinMM = mmhE
              , _hairpinLookup = hexloopE <> tetraloopE <> triloopE
              , _hairpinPenalty = endE
--              , _hairpinGGG = ddgDbg
--              , _hairpinCslope = ddgDbg
--              , _hairpinCintercept = ddgDbg
--              , _hairpinC3 = ddgDbg
              , _largeLoop = ddgDbg
              }
            intloopE = IntLoop
              { _intLoop1x1 = i11E
              , _intLoop1x2 = i21E
              , _intLoop2x2 = i22E
              , _intLoop1xn = mmi1nE
              , _intLoopMM = mmiE
              , _intLoopL = intE
              , _bulgeL = bulE
              }
            multiE = Multi
              { _mismatchMulti = mmmE
              }
            enthalpy = Turner2004
              { _hairpin = hairpinE
              , _stack   = sE
              , _intloop = intloopE
              , _multi   = multiE
              , _exterior = Exterior
                { _mismatchExterior = mmeE
                , _dangle5 = d5E
                , _dangle3 = d3E
                , _endPenalty = endE
                }
              }
          return (entropy,enthalpy)

-- | Parses a block of two stacked pairs.

blockPP :: String -> P (Dense Vector (Z:.R:.R:.R:.R) DDG)
blockPP s = string s *> someSpace *> (f <$> some ddg)
  where f = fromAssocs maxBound ddgDbg . zip pp

-- | Parser a block of a stacked pair, with two unpaired mismatched nucleotides.

blockPNN :: String -> P (Dense Vector (Z:.R:.R:.R:.R) DDG)
blockPNN s = string s *> someSpace *> (f <$> some ddg)
  where f = fromAssocs maxBound ddgDbg . zip pnn

blockPN :: String -> P (Dense Vector (Z:.R:.R:.R) DDG)
blockPN s = string s *> someSpace *> (f <$> some ddg)
  where f = fromAssocs maxBound ddgDbg . zip pn

blockPPNN :: String -> P (Dense Vector (Z:.R:.R:.R:.R:.R:.R) DDG)
blockPPNN s = string s *> someSpace *> (f <$> some ddg)
  where f = fromAssocs maxBound ddgDbg . zip ppnn

blockPPNNN :: String -> P (Dense Vector (Z:.R:.R:.R:.R:.R:.R:.R) DDG)
blockPPNNN s = string s *> someSpace *> (f <$> some ddg)
  where f = fromAssocs maxBound ddgDbg . zip ppnnn

blockPPNNNN :: String -> P (Dense Vector (Z:.R:.R:.R:.R:.R:.R:.R:.R) DDG)
blockPPNNNN s = string s *> someSpace *> (f <$> some ddg)
  where f = fromAssocs maxBound ddgDbg . zip ppnnnn

blockL :: String -> P (Dense Vector (Z:.Int) DDG)
blockL s = string s *> someSpace *> (f <$> some ddg)
  where f = fromAssocs (ZZ:..LtInt 31) ddgDbg . zip (map (Z:.) [0..30])

pmlLoop :: P (MlLoop DDG)
pmlLoop = MlLoop <$> ddg <*> ddg <*> ddg

loopMaps :: String -> P (M.Map (Vector R) DDG, M.Map (Vector R) DDG)
loopMaps s = do
  string s >> someSpace
  let f k v vE = (k,(v,vE))
      pkey = pNucV <* someSpace
  xs <- M.fromList <$> some (f <$> pkey <*> ddg <*> ddg)
  return (M.map fst xs, M.map snd xs)

ddg = DDG <$> fromInteger <$> (*100) <$> integer
    <|> (ddgDbg <$ string "INF" <* someSpace)

pNuc :: P R
pNuc = charRNA <$> TT.oneOf (map rnaChar acgu)

pNucV :: P (Vector R)
pNucV = fromList <$> some pNuc

pp = (\(a,b) (c,d) -> (Z:.a:.b:.c:.d)) <$> viennaPairsNN <*> viennaPairsNN
pnn = (\(a,b) c d  -> (Z:.a:.b:.c:.d)) <$> viennaPairsNN <*> nacgu <*> nacgu
pn = (\(a,b) c -> (Z:.a:.b:.c)) <$> viennaPairs <*> nacgu
ppnn = (\(a,b) (c,d) e f -> (Z:.a:.b:.c:.d:.e:.f)) <$> viennaPairsNN <*> viennaPairsNN <*> nacgu <*> nacgu
ppnnn = (\(a,b) (c,d) e f g -> (Z:.a:.b:.c:.d:.e:.f:.g)) <$> viennaPairsNN <*> viennaPairsNN <*> nacgu <*> nacgu <*> nacgu
ppnnnn = (\(a,b) (c,d) e f g h -> (Z:.a:.b:.c:.d:.e:.f:.g:.h)) <$> viennaPairs <*> viennaPairs <*> acgu <*> acgu <*> acgu <*> acgu

nacgu = N : acgu

ddgDbg = DDG 999999


-- * Simple test

runtest :: IO ()
runtest = do
  r <- runM $ runThrow $ fromFile "deps/BiobaseTurner/data/rna_turner2004.par"
  case r of
    Left (err::String) -> print err
    Right r -> print () >> print r

