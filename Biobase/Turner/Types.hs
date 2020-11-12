
-- | This module defines the @Turner 2004@ energy model data structure. The
-- model is composed of a large number of "lookup tables" which hold the
-- energies for different "loops".
--
-- TODO need a "Functor" instance over elements "e". Or alternatively, generic
-- programming to capture stuff going on in 'e'

module Biobase.Turner.Types where

import Control.Lens
import Data.Aeson (FromJSON, ToJSON)
import Data.Binary (Binary)
import Data.Data
import Data.Data.Lens
import Data.Default
import Data.Primitive.Types
import Data.Serialize (Serialize)
import Data.Typeable
import Data.Vector.Generic.Lens
import Data.Vector.Unboxed.Deriving
import Data.Vector.Unboxed (Vector)
import GHC.Generics
import qualified Data.Map as M
import qualified Data.Vector as V
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Unboxed as VU

import Biobase.Primary
import Biobase.Primary.Nuc.RNA
import Biobase.Secondary
import Biobase.Secondary.Basepair (Basepair)
import Biobase.Secondary.Vienna (ViennaPair)
import Biobase.Types.Energy
import Data.PrimitiveArray as PA
import qualified Biobase.Secondary.Basepair as SB
import qualified Biobase.Secondary.Vienna as SV



-- * Individual elements of the model.

-- | The @One@ loop or @Hairpin@. Should be monomeric.
--
-- @
-- 5'  ...PL...RQ...  3'
--        (.....)
-- @

data Hairpin vc ve c e = Hairpin
  { _hairpinLength      :: !(Dense ve (Z:.Int) e)
    -- ^ Contribution of the length of the unpaired region
  , _hairpinMM          :: !(Dense ve (Z:.c:.c:.c:.c) e)
    -- ^ The last match to first mismatch contribution. In 5'--3' order:
    -- @5' match - 5' mismatch - ... - 3' mismatch - 3' match@
  , _hairpinLookup      :: !(M.Map (vc c) e)
    -- ^ Tabulated energies for short hairpins that do not follow the generic
    -- scheme.
  , _hairpinAU          :: !(Dense ve (Z:.c:.c) e)
  , _largeLoop          :: !e
--  , _hairpinGGG         :: !e
--  , _hairpinCslope      :: !e
--  , _hairpinCintercept  :: !e
--  , _hairpinC3          :: !e
  }
  deriving (Generic)
makeLenses ''Hairpin

-- | Traversal over hairpin energies. This is polymorphic in the @ve@ vector-type for energies,
-- allowing for things like @show@ing the values.

hairpinE :: (VG.Vector ve a, VG.Vector we b) => Traversal (Hairpin vc ve c a) (Hairpin vc we c b) a b
{-# Inlinable hairpinE #-}
hairpinE f (Hairpin hl hmm lkup tAU llp) -- ggg sl ntr c3)
  = Hairpin <$> (denseV.vectorTraverse) f hl
            <*> (denseV.vectorTraverse) f hmm
            <*> traverse f lkup
            <*> (denseV.vectorTraverse) f tAU
            <*> f llp
--            <*> f ggg
--            <*> f sl
--            <*> f ntr
--            <*> f c3

deriving instance
  ( Show c, Show (LimitType c), Show e, Show (vc c), Show (ve e)
  ) => Show (Hairpin vc ve c e)

-- | Stacking helix contributions. Contains parameters for canonical stacks,
-- bulges, and interior loops.
--
-- TODO this is an @IntLoop@, too!
--
-- Order:
-- @
-- 5'   3'
-- GA...UC
-- 14   32
-- @
-- If drawn differently, one recognizes that pairs are ordered from 5' to 3'.
-- @
-- 3'    5'
--   A--U
--   G--C
-- 5'    3'
-- @

data Stack ve c e = Stack
  { _stacking :: !(Dense ve (Z:.c:.c:.c:.c) e)
  }

deriving instance (Show (ve e), Show c, Show e, Show (LimitType c)) => Show (Stack ve c e)

stackE :: (VG.Vector ve a, VG.Vector we b) => Traversal (Stack ve c a) (Stack we c b) a b
{-# Inlinable stackE #-}
stackE f (Stack s) = Stack <$> (denseV.vectorTraverse) f s

-- | All variants of interior loops and bulges.
--
-- Fix up the order of characters. Potentially, we want to introduce wrapping newtypes to prevent us
-- going crazy. Group by pairs?
--
-- TODO @_bulgeAU@ should be expanded to a full terminal mismatch system, even if the scoring is not
-- like this for now.

data IntLoop ve c e = IntLoop
  { _intLoop1x1 :: !(Dense ve (Z:.c:.c:.c:.c:.c:.c) e)
  , _intLoop1x2 :: !(Dense ve (Z:.c:.c:.c:.c:.c:.c:.c) e)
  , _intLoop2x2 :: !(Dense ve (Z:.c:.c:.c:.c:.c:.c:.c:.c) e)
  , _intLoop2x3 :: !(Dense ve (Z:.c:.c:.c:.c) e)
  , _intLoop1xn :: !(Dense ve (Z:.c:.c:.c:.c) e)
  , _intLoopMM  :: !(Dense ve (Z:.c:.c:.c:.c) e)
  , _intLoopL   :: !(Dense ve (Z:.Int) e)
  , _intLoopNinio :: !e
  , _intLoopMaxNinio :: !e
  , _bulgeL     :: !(Dense ve (Z:.Int) e)
  , _bulgeAU    :: !(Dense ve (Z:.c:.c) e)
  }
deriving instance (Show (ve e), Show c, Show e, Show (LimitType c)) => Show (IntLoop ve c e)

intLoopE :: (VG.Vector ve a, VG.Vector we b) => Traversal (IntLoop ve c a) (IntLoop we c b) a b
{-# Inlinable intLoopE #-}
intLoopE f (IntLoop i11 i12 i22 i23 i1n imm ilen ininio imaxninio blen bau)
  =   IntLoop
  <$> (denseV.vectorTraverse) f i11
  <*> (denseV.vectorTraverse) f i12
  <*> (denseV.vectorTraverse) f i22
  <*> (denseV.vectorTraverse) f i23
  <*> (denseV.vectorTraverse) f i1n
  <*> (denseV.vectorTraverse) f imm
  <*> (denseV.vectorTraverse) f ilen
  <*> f ininio
  <*> f imaxninio
  <*> (denseV.vectorTraverse) f blen
  <*> (denseV.vectorTraverse) f bau

-- * Multibranched parameters for loops.

data MlLoop e = MlLoop
  { _mlLoop :: !e
  , _mlOpen :: !e
  , _mlNuc  :: !e
  }

-- * The full model.

data Turner2004 vc ve c e = Turner2004
  { _hairpin :: !(Hairpin vc ve c e)
    -- ^ Hairpin contributions
  , _stack   :: !(Stack ve c e)
  }
makeLenses ''Turner2004

deriving instance (Show c, Show e, Show (vc c), Show (ve e), Show (LimitType c)) => Show (Turner2004 vc ve c e)

-- | Traverse all energy fields with a polymorphic and energy-vector changing traversal.

turner2004E :: (VG.Vector ve a, VG.Vector we b) => Traversal (Turner2004 vc ve c a) (Turner2004 vc we c b) a b
{-# Inlinable turner2004E #-}
turner2004E f (Turner2004 hp s) = Turner2004 <$> hairpinE f hp <*> stackE f s

