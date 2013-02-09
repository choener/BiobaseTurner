{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE TypeOperators #-}

-- | The 'Turner2004' data structure reflects the RNA (and DNA) energy
-- parameters known as the Turner 2004 data set.
--
-- In general, have a look here:
-- <http://rna.urmc.rochester.edu/NNDB/turner04/index.html> where parameters
-- are explained.
--
-- TODO need a "Functor" instance over elements "e". Or alternatively, generic
-- programming to capture stuff going on in 'e'

module Biobase.Turner where


import Control.Lens
import Data.Array.Repa.Index
import Data.ByteString
import Data.Map as M
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Generic.Mutable as VGM
import Data.Primitive.Types

import Biobase.Primary
import Biobase.Secondary
import Data.PrimitiveArray
import Data.PrimitiveArray.Zero



-- | The actual Turner parameters return energies in Double format.

newtype Energy = Energy Double
  deriving (Eq,Ord,Num,Read,Show)

deriving instance Prim Energy
deriving instance VGM.MVector VU.MVector Energy
deriving instance VG.Vector   VU.Vector  Energy
deriving instance VU.Unbox Energy

-- | The Turner model with 'Energy's.

type Turner2004 = Turner2004Model VU.Vector Energy

-- | The Turner energy tables. Parametrized over the storing vector type 'v'
-- and the actual element type 'e'.

data Turner2004Model v e = Turner2004Model
  { _stack              :: A PP v e
  , _dangle3            :: A PN v e
  , _dangle5            :: A PN v e
  , _hairpinL           :: A DIM1 v e
  , _hairpinMM          :: A PNN v e
  , _hairpinLookup      :: M.Map ByteString e
  , _hairpinGGG         :: e
  , _hairpinCslope      :: e
  , _hairpinCintercept  :: e
  , _hairpinC3          :: e
  , _bulgeL             :: A DIM1 v e
  , _bulgeSingleC       :: e
  , _iloop1x1           :: A PPNN v e
  , _iloop2x1           :: A PPNNN v e
  , _iloop2x2           :: A PPNNNN v e
  , _iloopMM            :: A PNN v e
  , _iloop2x3MM         :: A PNN v e
  , _iloop1xnMM         :: A PNN v e
  , _iloopL             :: A DIM1 v e
  , _multiMM            :: A PNN v e
  , _ninio              :: e
  , _maxNinio           :: e
  , _multiOffset        :: e
  , _multiNuc           :: e
  , _multiHelix         :: e
  , _multiAsym          :: e
  , _multiStrain        :: e
  , _extMM              :: A PNN v e
  , _coaxial            :: A PP v e -- no intervening unpaired nucleotides
  , _coaxStack          :: A PNN v e
  , _tStackCoax         :: A PNN v e
  , _largeLoop          :: e
  , _termAU             :: e
  , _intermolecularInit :: e
  } deriving (Show)

type PP = (Z:.Nuc:.Nuc:.Nuc:.Nuc)
type PN = (Z:.Nuc:.Nuc:.Nuc)
type PNN = (Z:.Nuc:.Nuc:.Nuc:.Nuc)
type PPNN = PP:.Nuc:.Nuc
type PPNNN = PPNN:.Nuc
type PPNNNN = PPNNN:.Nuc

makeLenses ''Turner2004Model

