
-- | This module defines the @Turner 2004@ energy model data structure. The
-- model is composed of a large number of "lookup tables" which hold the
-- energies for different "loops".
--
-- TODO need a "Functor" instance over elements "e". Or alternatively, generic
-- programming to capture stuff going on in 'e'

module Biobase.Turner.Types where

import           Control.Lens
import           Data.Primitive.Types
import           Data.Vector.Unboxed.Deriving
import           Data.Vector.Unboxed (Vector)
import qualified Data.Map as M
import qualified Data.Vector.Unboxed as VU
import           GHC.Generics
import           Data.Binary (Binary)
import           Data.Serialize (Serialize)
import           Data.Aeson (FromJSON, ToJSON)
import           Data.Default

import           Biobase.Primary
import           Biobase.Primary.Nuc.RNA
import           Biobase.Secondary
import           Biobase.Secondary.Basepair (Basepair)
import           Biobase.Secondary.Vienna (ViennaPair)
import           Biobase.Types.Energy
import           Data.PrimitiveArray as PA
import qualified Biobase.Secondary.Basepair as SB
import qualified Biobase.Secondary.Vienna as SV



-- | The Turner model with 'Energy's.

type Turner2004 = Turner2004Model Basepair DG

type Vienna2004 = Turner2004Model ViennaPair DDG

type PU p u = Z:.p:.u
type PP p   = Z:.p:.p
type PUU p u = PU p u :.u
type PPUU p u = PP p:.u:.u
type PPU3 p u = PPUU p u:.u
type PPU4 p u = PPU3 p u:.u

-- | The Turner energy tables. Parametrized over the actual element type
-- 'e'.

data Turner2004Model p u e = Turner2004Model
  { _stack, _coaxial    ∷ !(Unboxed (PP p) e)
  , _dangle3, _dangle5  ∷ !(Unboxed (PU p u) e)
  , _hairpinL           ∷ !(Vector e)
  , _hairpinMM          ∷ !(Unboxed (Z:.p:.u:.u) e)
  , _hairpinLookup      ∷ !(M.Map (Primary u) e)
  , _hairpinGGG         ∷ !e
  , _hairpinCslope      ∷ !e
  , _hairpinCintercept  ∷ !e
  , _hairpinC3          ∷ !e
  , _bulgeL             ∷ !(Vector e)
  , _bulgeSingleC       ∷ !e
  , _iloop1x1           ∷ !(Unboxed (PPUU p u) e)
  , _iloop2x1           ∷ !(Unboxed (PPU3 p u) e)
  , _iloop2x2           ∷ !(Unboxed (PPU4 p u) e)
  , _iloopL             ∷ !(Vector e)
  , _iloopMM, _iloop2x3MM, _iloop1xnMM, _multiMM, _exteriorMM, _coaxStack, _tStackCoax ∷ !(Unboxed (PUU p u) e)
  , _ninio, _maxNinio, _multiOffset, _multiNuc, _multiHelix, _multiAsym, _multiStrain, _largeLoop, _termAU, _intermolecularInit ∷ !e
  } deriving (Generic)

type N = Letter RNA

{-
type PP = (Z:.Letter RNA:.Letter RNA:.Letter RNA:.Letter RNA)
type PN = (Z:.Letter RNA:.Letter RNA:.Letter RNA)
type PNN = (Z:.Letter RNA:.Letter RNA:.Letter RNA:.Letter RNA)
type PPNN = PP:.Letter RNA:.Letter RNA
type PPNNN = PPNN:.Letter RNA
type PPNNNN = PPNNN:.Letter RNA
-}
type VPP = (Z:.ViennaPair:.ViennaPair)
type VPN = (Z:.ViennaPair:.Letter RNA)
type VPNN = VPN:.Letter RNA
type VPPNN = VPP:.Letter RNA:.Letter RNA
type VPPNNN = VPPNN:.Letter RNA
type VPPNNNN = VPPNNN:.Letter RNA

makeLenses ''Turner2004Model

--instance (Binary p, VU.Unbox e, Binary e                                        ) => Binary    (Turner2004Model p e)
--instance (Serialize p, VU.Unbox e, Serialize e                                  ) => Serialize (Turner2004Model p e)
--instance (FromJSON p, VU.Unbox e, FromJSON e , FromJSON (M.Map (Primary RNA) e) ) => FromJSON  (Turner2004Model p e)
--instance (ToJSON p, VU.Unbox e, ToJSON e   , ToJSON   (M.Map (Primary RNA) e)   ) => ToJSON    (Turner2004Model p e)



-- | Map a function over all 'e' elements.

emap
  ∷ (PA.Index p, PA.Index u, VU.Unbox e, VU.Unbox e')
  ⇒ (e → e')
  -- ^ conversion of "energies" @e@ to @e'@ with potentially different type.
  -- the @e@'s could very well be probabilities.
  → Turner2004Model p u e
  -- the input turner model. We may not change the paired @p@ or unpaired @u@
  -- types, only @e@ with 'emap'.
  → Turner2004Model p u e'
emap f Turner2004Model{..} = Turner2004Model
  { _stack              = PA.map f _stack
  , _dangle3            = PA.map f _dangle3
  , _dangle5            = PA.map f _dangle5
  , _hairpinL           = VU.map f _hairpinL
  , _hairpinMM          = PA.map f _hairpinMM
  , _hairpinLookup      = M.map f _hairpinLookup
  , _hairpinGGG         = f _hairpinGGG
  , _hairpinCslope      = f _hairpinCslope
  , _hairpinCintercept  = f _hairpinCintercept
  , _hairpinC3          = f _hairpinC3
  , _bulgeL             = VU.map f _bulgeL
  , _bulgeSingleC       = f _bulgeSingleC
  , _iloop1x1           = PA.map f _iloop1x1
  , _iloop2x1           = PA.map f _iloop2x1
  , _iloop2x2           = PA.map f _iloop2x2
  , _iloopMM            = PA.map f _iloopMM
  , _iloop2x3MM         = PA.map f _iloop2x3MM
  , _iloop1xnMM         = PA.map f _iloop1xnMM
  , _iloopL             = VU.map f _iloopL
  , _multiMM            = PA.map f _multiMM
  , _ninio              = f _ninio
  , _maxNinio           = f _maxNinio
  , _multiOffset        = f _multiOffset
  , _multiNuc           = f _multiNuc
  , _multiHelix         = f _multiHelix
  , _multiAsym          = f _multiAsym
  , _multiStrain        = f _multiStrain
  , _exteriorMM              = PA.map f _exteriorMM
  , _coaxial            = PA.map f _coaxial
  , _coaxStack          = PA.map f _coaxStack
  , _tStackCoax         = PA.map f _tStackCoax
  , _largeLoop          = f _largeLoop
  , _termAU             = f _termAU
  , _intermolecularInit = f _intermolecularInit
  }

-- | An empty model

{-
emptyTurnerModel :: (Default e, VU.Unbox e) => Turner2004Model Basepair e
emptyTurnerModel = Turner2004Model
  { _stack              = fromAssocs minBPP def []
  , _dangle3            = fromAssocs minBPB def []
  , _dangle5            = fromAssocs minBPB def []
  , _hairpinL           = VU.empty
  , _hairpinMM          = fromAssocs minBPBB def []
  , _hairpinLookup      = M.empty
  , _hairpinGGG         = def
  , _hairpinCslope      = def
  , _hairpinCintercept  = def
  , _hairpinC3          = def
  , _bulgeL             = VU.empty
  , _bulgeSingleC       = def
  , _iloop1x1           = fromAssocs minBPPBB def []
  , _iloop2x1           = fromAssocs minBPPBBB def []
  , _iloop2x2           = fromAssocs minBPPBBBB def []
  , _iloopMM            = fromAssocs minBPBB def []
  , _iloop2x3MM         = fromAssocs minBPBB def []
  , _iloop1xnMM         = fromAssocs minBPBB def []
  , _iloopL             = VU.empty
  , _multiMM            = fromAssocs minBPBB def []
  , _ninio              = def
  , _maxNinio           = def
  , _multiOffset        = def
  , _multiNuc           = def
  , _multiHelix         = def
  , _multiAsym          = def
  , _multiStrain        = def
  , _exteriorMM         = fromAssocs minBPBB def []
  , _coaxial            = fromAssocs minBPP def []
  , _coaxStack          = fromAssocs minBPBB def []
  , _tStackCoax         = fromAssocs minBPBB def []
  , _largeLoop          = def
  , _termAU             = def
  , _intermolecularInit = def
  }

emptyViennaModel :: (Default e, VU.Unbox e) => Turner2004Model ViennaPair e
emptyViennaModel = Turner2004Model
  { _stack              = fromAssocs minVPP def []
  , _dangle3            = fromAssocs minVPB def []
  , _dangle5            = fromAssocs minVPB def []
  , _hairpinL           = VU.empty
  , _hairpinMM          = fromAssocs minVPBB def []
  , _hairpinLookup      = M.empty
  , _hairpinGGG         = def
  , _hairpinCslope      = def
  , _hairpinCintercept  = def
  , _hairpinC3          = def
  , _bulgeL             = VU.empty
  , _bulgeSingleC       = def
  , _iloop1x1           = fromAssocs minVPPBB def []
  , _iloop2x1           = fromAssocs minVPPBBB def []
  , _iloop2x2           = fromAssocs minVPPBBBB def []
  , _iloopMM            = fromAssocs minVPBB def []
  , _iloop2x3MM         = fromAssocs minVPBB def []
  , _iloop1xnMM         = fromAssocs minVPBB def []
  , _iloopL             = VU.empty
  , _multiMM            = fromAssocs minVPBB def []
  , _ninio              = def
  , _maxNinio           = def
  , _multiOffset        = def
  , _multiNuc           = def
  , _multiHelix         = def
  , _multiAsym          = def
  , _multiStrain        = def
  , _exteriorMM         = fromAssocs minVPBB def []
  , _coaxial            = fromAssocs minVPP def []
  , _coaxStack          = fromAssocs minVPBB def []
  , _tStackCoax         = fromAssocs minVPBB def []
  , _largeLoop          = def
  , _termAU             = def
  , _intermolecularInit = def
  }
-}

minBPP     = Z:.SB.AA:.SB.AA
maxBPP     = Z:.SB.NS:.SB.NS
minBP      = Z:.SB.AA
maxBP      = Z:.SB.NS
minBPB     = minBP:.A
maxBPB     = maxBP:.N
minBPBB    = minBPB:.A
maxBPBB    = maxBPB:.N
minBPPBB   = minBPP:.A:.A
maxBPPBB   = maxBPP:.N:.N
minBPPBBB  = minBPPBB:.A
maxBPPBBB  = maxBPPBB:.N
minBPPBBBB = minBPPBBB:.A
maxBPPBBBB = maxBPPBBB:.N

minVPP     = Z:.SV.NP:.SV.NP
maxVPP     = Z:.SV.NS:.SV.NS
minVP      = Z:.SV.NP
maxVP      = Z:.SV.NS
minVPB     = minVP:.A
maxVPB     = maxVP:.N
minVPBB    = minVPB:.A
maxVPBB    = maxVPB:.N
minVPPBB   = minVPP:.A:.A
maxVPPBB   = maxVPP:.N:.N
minVPPBBB  = minVPPBB:.A
maxVPPBBB  = maxVPPBB:.N
minVPPBBBB = minVPPBBB:.A
maxVPPBBBB = maxVPPBBB:.N

