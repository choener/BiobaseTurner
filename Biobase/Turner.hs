{-# LANGUAGE TypeOperators #-}

-- | The 'Turner2004' data structure reflects the RNA (and DNA) energy
-- parameters known as the Turner 2004 data set.
--
-- In general, have a look here:
-- <http://rna.urmc.rochester.edu/NNDB/turner04/index.html> where parameters
-- are explained.

module Biobase.Turner where

import Data.ByteString
import Data.Map as M
import Data.Array.Repa.Index

import Biobase.Primary
import Biobase.Secondary
import Data.PrimitiveArray
import Data.PrimitiveArray.Unboxed

-- | The parameters. Turner parameters are set by the Import module for
-- nucleotides n,a,c,g,u. All values that are not read (or are ".") will end up
-- with a value > 100K.
--
-- TODO use 'Energy' instead of 'Double'
--
-- TODO specialized shape types for pairs?

type PP = (Z:.Nuc:.Nuc:.Nuc:.Nuc)
type PN = (Z:.Nuc:.Nuc:.Nuc)
type PNN = (Z:.Nuc:.Nuc:.Nuc:.Nuc)
type PPNN = PP:.Nuc:.Nuc
type PPNNN = PPNN:.Nuc
type PPNNNN = PPNNN:.Nuc

data Turner2004 = Turner2004
  { stack :: PrimArray PP Double
  , dangle3 :: PrimArray PN Double
  , dangle5 :: PrimArray PN Double
  , hairpinL :: PrimArray DIM1 Double
  , hairpinMM :: PrimArray PNN Double
  , hairpinLookup :: M.Map ByteString Double
  , hairpinGGG :: Double
  , hairpinCslope :: Double
  , hairpinCintercept :: Double
  , hairpinC3 :: Double
  , bulgeL :: PrimArray DIM1 Double
  , bulgeSingleC :: Double
  , iloop1x1 :: PrimArray PPNN Double
  , iloop2x1 :: PrimArray PPNNN Double
  , iloop2x2 :: PrimArray PPNNNN Double
  , iloopMM :: PrimArray PNN Double
  , iloop2x3MM :: PrimArray PNN Double
  , iloop1xnMM :: PrimArray PNN Double
  , iloopL :: PrimArray DIM1 Double
  , multiMM :: PrimArray PNN Double
  , ninio :: Double
  , maxNinio :: Double
  , multiOffset :: Double
  , multiNuc :: Double
  , multiHelix :: Double
  , multiAsym :: Double
  , multiStrain :: Double
  , extMM :: PrimArray PNN Double
  , coaxial :: PrimArray PP Double -- no intervening unpaired nucleotides
  , coaxStack :: PrimArray PNN Double
  , tStackCoax :: PrimArray PNN Double
  , largeLoop :: Double
  , termAU :: Double
  , intermolecularInit :: Double
  } deriving (Show)

