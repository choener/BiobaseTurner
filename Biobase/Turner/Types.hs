
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

type Turner2004 = Turner2004Model Basepair DeltaGibbs

type Vienna2004 = Turner2004Model ViennaPair DeltaDekaGibbs

-- | The Turner energy tables. Parametrized over the actual element type
-- 'e'.

data Turner2004Model p e = Turner2004Model
  { _stack              :: !(Unboxed (Z:.p:.p) e)
  , _dangle3            :: !(Unboxed (Z:.p:.N) e)
  , _dangle5            :: !(Unboxed (Z:.p:.N) e)
  , _hairpinL           :: !(Vector e)
  , _hairpinMM          :: !(Unboxed (Z:.p:.N:.N) e)
  , _hairpinLookup      :: !(M.Map (Primary RNA) e)
  , _hairpinGGG         :: !e
  , _hairpinCslope      :: !e
  , _hairpinCintercept  :: !e
  , _hairpinC3          :: !e
  , _bulgeL             :: !(Vector e)
  , _bulgeSingleC       :: !e
  , _iloop1x1           :: !(Unboxed (Z:.p:.p:.N:.N) e)
  , _iloop2x1           :: !(Unboxed (Z:.p:.p:.N:.N:.N) e)
  , _iloop2x2           :: !(Unboxed (Z:.p:.p:.N:.N:.N:.N) e)
  , _iloopMM            :: !(Unboxed (Z:.p:.N:.N) e)
  , _iloop2x3MM         :: !(Unboxed (Z:.p:.N:.N) e)
  , _iloop1xnMM         :: !(Unboxed (Z:.p:.N:.N) e)
  , _iloopL             :: !(Vector e)
  , _multiMM            :: !(Unboxed (Z:.p:.N:.N) e)
  , _ninio              :: !e
  , _maxNinio           :: !e
  , _multiOffset        :: !e
  , _multiNuc           :: !e
  , _multiHelix         :: !e
  , _multiAsym          :: !e
  , _multiStrain        :: !e
  , _exteriorMM         :: !(Unboxed (Z:.p:.N:.N) e)
  , _coaxial            :: !(Unboxed (Z:.p:.p) e) -- no intervening unpaired nucleotides
  , _coaxStack          :: !(Unboxed (Z:.p:.N:.N) e)
  , _tStackCoax         :: !(Unboxed (Z:.p:.N:.N) e)
  , _largeLoop          :: !e
  , _termAU             :: !e
  , _intermolecularInit :: !e
  } deriving (Show,Generic)

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

instance (Binary p, VU.Unbox e, Binary e                                        ) => Binary    (Turner2004Model p e)
instance (Serialize p, VU.Unbox e, Serialize e                                  ) => Serialize (Turner2004Model p e)
instance (FromJSON p, VU.Unbox e, FromJSON e , FromJSON (M.Map (Primary RNA) e) ) => FromJSON  (Turner2004Model p e)
instance (ToJSON p, VU.Unbox e, ToJSON e   , ToJSON   (M.Map (Primary RNA) e)   ) => ToJSON    (Turner2004Model p e)



-- | Map a function over all 'e' elements.

emap :: (PA.Index p, VU.Unbox e, VU.Unbox e') => (e -> e') -> Turner2004Model p e -> Turner2004Model p e'
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

emptyTurnerModel :: (Default e, VU.Unbox e) => Turner2004Model Basepair e
emptyTurnerModel = Turner2004Model
  { _stack              = fromAssocs minBPP minBPP def []
  , _dangle3            = fromAssocs minBPB minBPB def []
  , _dangle5            = fromAssocs minBPB minBPB def []
  , _hairpinL           = VU.empty
  , _hairpinMM          = fromAssocs minBPBB minBPBB def []
  , _hairpinLookup      = M.empty
  , _hairpinGGG         = def
  , _hairpinCslope      = def
  , _hairpinCintercept  = def
  , _hairpinC3          = def
  , _bulgeL             = VU.empty
  , _bulgeSingleC       = def
  , _iloop1x1           = fromAssocs minBPPBB minBPPBB def []
  , _iloop2x1           = fromAssocs minBPPBBB minBPPBBB def []
  , _iloop2x2           = fromAssocs minBPPBBBB minBPPBBBB def []
  , _iloopMM            = fromAssocs minBPBB minBPBB def []
  , _iloop2x3MM         = fromAssocs minBPBB minBPBB def []
  , _iloop1xnMM         = fromAssocs minBPBB minBPBB def []
  , _iloopL             = VU.empty
  , _multiMM            = fromAssocs minBPBB minBPBB def []
  , _ninio              = def
  , _maxNinio           = def
  , _multiOffset        = def
  , _multiNuc           = def
  , _multiHelix         = def
  , _multiAsym          = def
  , _multiStrain        = def
  , _exteriorMM         = fromAssocs minBPBB minBPBB def []
  , _coaxial            = fromAssocs minBPP minBPP def []
  , _coaxStack          = fromAssocs minBPBB minBPBB def []
  , _tStackCoax         = fromAssocs minBPBB minBPBB def []
  , _largeLoop          = def
  , _termAU             = def
  , _intermolecularInit = def
  }

emptyViennaModel :: (Default e, VU.Unbox e) => Turner2004Model ViennaPair e
emptyViennaModel = Turner2004Model
  { _stack              = fromAssocs minVPP minVPP def []
  , _dangle3            = fromAssocs minVPB minVPB def []
  , _dangle5            = fromAssocs minVPB minVPB def []
  , _hairpinL           = VU.empty
  , _hairpinMM          = fromAssocs minVPBB minVPBB def []
  , _hairpinLookup      = M.empty
  , _hairpinGGG         = def
  , _hairpinCslope      = def
  , _hairpinCintercept  = def
  , _hairpinC3          = def
  , _bulgeL             = VU.empty
  , _bulgeSingleC       = def
  , _iloop1x1           = fromAssocs minVPPBB minVPPBB def []
  , _iloop2x1           = fromAssocs minVPPBBB minVPPBBB def []
  , _iloop2x2           = fromAssocs minVPPBBBB minVPPBBBB def []
  , _iloopMM            = fromAssocs minVPBB minVPBB def []
  , _iloop2x3MM         = fromAssocs minVPBB minVPBB def []
  , _iloop1xnMM         = fromAssocs minVPBB minVPBB def []
  , _iloopL             = VU.empty
  , _multiMM            = fromAssocs minVPBB minVPBB def []
  , _ninio              = def
  , _maxNinio           = def
  , _multiOffset        = def
  , _multiNuc           = def
  , _multiHelix         = def
  , _multiAsym          = def
  , _multiStrain        = def
  , _exteriorMM         = fromAssocs minVPBB minVPBB def []
  , _coaxial            = fromAssocs minVPP minVPP def []
  , _coaxStack          = fromAssocs minVPBB minVPBB def []
  , _tStackCoax         = fromAssocs minVPBB minVPBB def []
  , _largeLoop          = def
  , _termAU             = def
  , _intermolecularInit = def
  }

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

