
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
import           Biobase.Types.Energy
import           Data.PrimitiveArray as PA



-- | The Turner model with 'Energy's.

type Turner2004 = Turner2004Model DeltaGibbs

type Vienna2004 = Turner2004Model DeltaDekaGibbs

-- | The Turner energy tables. Parametrized over the actual element type
-- 'e'.

data Turner2004Model e = Turner2004Model
  { _stack              :: !(Unboxed PP e)
  , _dangle3            :: !(Unboxed PN e)
  , _dangle5            :: !(Unboxed PN e)
  , _hairpinL           :: !(Vector e)
  , _hairpinMM          :: !(Unboxed PNN e)
  , _hairpinLookup      :: !(M.Map (Primary RNA) e)
  , _hairpinGGG         :: !e
  , _hairpinCslope      :: !e
  , _hairpinCintercept  :: !e
  , _hairpinC3          :: !e
  , _bulgeL             :: !(Vector e)
  , _bulgeSingleC       :: !e
  , _iloop1x1           :: !(Unboxed PPNN e)
  , _iloop2x1           :: !(Unboxed PPNNN e)
  , _iloop2x2           :: !(Unboxed PPNNNN e)
  , _iloopMM            :: !(Unboxed PNN e)
  , _iloop2x3MM         :: !(Unboxed PNN e)
  , _iloop1xnMM         :: !(Unboxed PNN e)
  , _iloopL             :: !(Vector e)
  , _multiMM            :: !(Unboxed PNN e)
  , _ninio              :: !e
  , _maxNinio           :: !e
  , _multiOffset        :: !e
  , _multiNuc           :: !e
  , _multiHelix         :: !e
  , _multiAsym          :: !e
  , _multiStrain        :: !e
  , _exteriorMM         :: !(Unboxed PNN e)
  , _coaxial            :: !(Unboxed PP e) -- no intervening unpaired nucleotides
  , _coaxStack          :: !(Unboxed PNN e)
  , _tStackCoax         :: !(Unboxed PNN e)
  , _largeLoop          :: !e
  , _termAU             :: !e
  , _intermolecularInit :: !e
  } deriving (Show,Generic)

type PP = (Z:.Letter RNA:.Letter RNA:.Letter RNA:.Letter RNA)
type PN = (Z:.Letter RNA:.Letter RNA:.Letter RNA)
type PNN = (Z:.Letter RNA:.Letter RNA:.Letter RNA:.Letter RNA)
type PPNN = PP:.Letter RNA:.Letter RNA
type PPNNN = PPNN:.Letter RNA
type PPNNNN = PPNNN:.Letter RNA

makeLenses ''Turner2004Model

instance (VU.Unbox e, Binary e                                      ) => Binary    (Turner2004Model e)
instance (VU.Unbox e, Serialize e                                   ) => Serialize (Turner2004Model e)
instance (VU.Unbox e, FromJSON e , FromJSON (M.Map (Primary RNA) e) ) => FromJSON  (Turner2004Model e)
instance (VU.Unbox e, ToJSON e   , ToJSON   (M.Map (Primary RNA) e) ) => ToJSON    (Turner2004Model e)



-- | Map a function over all 'e' elements.

emap :: (VU.Unbox e, VU.Unbox e') => (e -> e') -> Turner2004Model e -> Turner2004Model e'
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

emptyModel :: (Default e, VU.Unbox e) => Turner2004Model e
emptyModel = Turner2004Model
  { _stack              = fromAssocs minPP minPP def []
  , _dangle3            = fromAssocs minPB minPB def []
  , _dangle5            = fromAssocs minPB minPB def []
  , _hairpinL           = VU.empty
  , _hairpinMM          = fromAssocs minPBB minPBB def []
  , _hairpinLookup      = M.empty
  , _hairpinGGG         = def
  , _hairpinCslope      = def
  , _hairpinCintercept  = def
  , _hairpinC3          = def
  , _bulgeL             = VU.empty
  , _bulgeSingleC       = def
  , _iloop1x1           = fromAssocs minPPBB minPPBB def []
  , _iloop2x1           = fromAssocs minPPBBB minPPBBB def []
  , _iloop2x2           = fromAssocs minPPBBBB minPPBBBB def []
  , _iloopMM            = fromAssocs minPBB minPBB def []
  , _iloop2x3MM         = fromAssocs minPBB minPBB def []
  , _iloop1xnMM         = fromAssocs minPBB minPBB def []
  , _iloopL             = VU.empty
  , _multiMM            = fromAssocs minPBB minPBB def []
  , _ninio              = def
  , _maxNinio           = def
  , _multiOffset        = def
  , _multiNuc           = def
  , _multiHelix         = def
  , _multiAsym          = def
  , _multiStrain        = def
  , _exteriorMM         = fromAssocs minPBB minPBB def []
  , _coaxial            = fromAssocs minPP minPP def []
  , _coaxStack          = fromAssocs minPBB minPBB def []
  , _tStackCoax         = fromAssocs minPBB minPBB def []
  , _largeLoop          = def
  , _termAU             = def
  , _intermolecularInit = def
  }

minPP     = Z:.N:.N:.N:.N -- (minP,minP)
maxPP     = Z:.U:.U:.U:.U -- (maxP,maxP)
minP      = Z:.N:.N -- (nN,nN)
maxP      = Z:.U:.U -- (nU,nU)
minPB     = minP:.N -- (minP,nN)
maxPB     = maxP:.U -- (maxP,nU)
minPBB    = minPB:.N -- (minP,nN,nN)
maxPBB    = maxPB:.U -- (maxP,nU,nU)
minPPBB   = minPP:.N:.N -- (minP,minP,(nN,nN))
maxPPBB   = maxPP:.U:.U -- (maxP,maxP,(nU,nU))
minPPBBB  = minPPBB:.N -- (minP,minP,(nN,nN,nN))
maxPPBBB  = maxPPBB:.U -- (maxP,maxP,(nU,nU,nU))
minPPBBBB = minPPBBB:.N -- (minP,minP,(nN,nN,nN,nN))
maxPPBBBB = maxPPBBB:.U -- (maxP,maxP,(nU,nU,nU,nU))

