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
import Data.PrimitiveArray.Unboxed.Zero

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
  { stack :: Arr0 PP Double
  , dangle3 :: Arr0 PN Double
  , dangle5 :: Arr0 PN Double
  , hairpinL :: Arr0 DIM1 Double
  , hairpinMM :: Arr0 PNN Double
  , hairpinLookup :: M.Map ByteString Double
  , hairpinGGG :: Double
  , hairpinCslope :: Double
  , hairpinCintercept :: Double
  , hairpinC3 :: Double
  , bulgeL :: Arr0 DIM1 Double
  , bulgeSingleC :: Double
  , iloop1x1 :: Arr0 PPNN Double
  , iloop2x1 :: Arr0 PPNNN Double
  , iloop2x2 :: Arr0 PPNNNN Double
  , iloopMM :: Arr0 PNN Double
  , iloop2x3MM :: Arr0 PNN Double
  , iloop1xnMM :: Arr0 PNN Double
  , iloopL :: Arr0 DIM1 Double
  , multiMM :: Arr0 PNN Double
  , ninio :: Double
  , maxNinio :: Double
  , multiOffset :: Double
  , multiNuc :: Double
  , multiHelix :: Double
  , multiAsym :: Double
  , multiStrain :: Double
  , extMM :: Arr0 PNN Double
  , coaxial :: Arr0 PP Double -- no intervening unpaired nucleotides
  , coaxStack :: Arr0 PNN Double
  , tStackCoax :: Arr0 PNN Double
  , largeLoop :: Double
  , termAU :: Double
  , intermolecularInit :: Double
  } deriving ()

