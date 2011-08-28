
-- | The 'Turner2004' data structure reflects the RNA (and DNA) energy
-- parameters known as the Turner 2004 data set.
--
-- In general, have a look here:
-- <http://rna.urmc.rochester.edu/NNDB/turner04/index.html> where parameters
-- are explained.

module Biobase.Turner where

import Data.ByteString
import Data.Map as M

import Biobase.Primary
import Biobase.Secondary
import Data.PrimitiveArray
import Data.PrimitiveArray.Ix

-- | The parameters. Turner parameters are set by the Import module for
-- nucleotides n,a,c,g,u. All values that are not read (or are ".") will end up
-- with a value > 100K.
--
-- [1] Yes, such instances are easily created, but I don't want to pull in
-- another library and I don't want to create an Ix instance here.
--
-- TODO use 'Energy' instead of 'Double'

data Turner2004 = Turner2004
  { stack :: PrimArray (Pair,Pair) Double
  , dangle3 :: PrimArray PN Double
  , dangle5 :: PrimArray PN Double
  , hairpinL :: PrimArray Int Double
  , hairpinMM :: PrimArray PNN Double
  , hairpinLookup :: M.Map ByteString Double
  , hairpinGGG :: Double
  , hairpinCslope :: Double
  , hairpinCintercept :: Double
  , hairpinC3 :: Double
  , bulgeL :: PrimArray Int Double
  , bulgeSingleC :: Double
  , iloop1x1 :: PrimArray (Pair,Pair,(Nuc,Nuc)) Double
  , iloop2x1 :: PrimArray (Pair,Pair,(Nuc,Nuc,Nuc)) Double
  , iloop2x2 :: PrimArray (Pair,Pair,(Nuc,Nuc,Nuc,Nuc)) Double -- yeah, 6-tuple ix instances :-( [1]
  , iloopMM :: PrimArray PNN Double
  , iloop2x3MM :: PrimArray PNN Double
  , iloop1xnMM :: PrimArray PNN Double
  , iloopL :: PrimArray Int Double
  , multiMM :: PrimArray PNN Double
  , ninio :: Double
  , maxNinio :: Double
  , multiOffset :: Double
  , multiNuc :: Double
  , multiHelix :: Double
  , multiAsym :: Double
  , multiStrain :: Double
  , extMM :: PrimArray PNN Double
  , coaxial :: PrimArray (Pair,Pair) Double -- no intervening unpaired nucleotides
  , coaxStack :: PrimArray PNN Double
  , tStackCoax :: PrimArray PNN Double
  , largeLoop :: Double
  , termAU :: Double
  , intermolecularInit :: Double
  } deriving (Read,Show)

type PNN = (Pair,Nuc,Nuc)
type PN  = (Pair,Nuc)
