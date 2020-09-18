
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
import           Data.Typeable
import           Data.Data
import           Data.Data.Lens
import           Data.Vector.Generic.Lens
import qualified Data.Vector.Generic as VG
import qualified Data.Vector as V

import           Biobase.Primary
import           Biobase.Primary.Nuc.RNA
import           Biobase.Secondary
import           Biobase.Secondary.Basepair (Basepair)
import           Biobase.Secondary.Vienna (ViennaPair)
import           Biobase.Types.Energy
import           Data.PrimitiveArray as PA
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
  { _hairpinLength      :: !(ve e)
    -- ^ Contribution of the length of the unpaired region
  , _hairpinMM          :: !(Dense ve (Z:.c:.c:.c:.c) e)
    -- ^ The last match to first mismatch contribution. In 5'--3' order:
    -- @5' match - 5' mismatch - ... - 3' mismatch - 3' match@
  , _hairpinLookup      :: !(M.Map (vc c) e)
    -- ^ Tabulated energies for short hairpins that do not follow the generic
    -- scheme.
  , _hairpinGGG         :: !e
  , _hairpinCslope      :: !e
  , _hairpinCintercept  :: !e
  , _hairpinC3          :: !e
  , _hairpinTerm        :: !(Dense ve (Z:.c:.c) e)
  }
  deriving (Generic)
makeLenses ''Hairpin

-- | Traversal over hairpin energies. This is polymorphic in the @ve@ vector-type for energies,
-- allowing for things like @show@ing the values.

hairpinE :: (VG.Vector ve a, VG.Vector we b) => Traversal (Hairpin vc ve c a) (Hairpin vc we c b) a b
{-# Inlinable hairpinE #-}
hairpinE f (Hairpin hl hmm lkup ggg sl ntr c3 ht)
  = Hairpin <$> vectorTraverse f hl
            <*> (denseV.vectorTraverse) f hmm
            <*> traverse f lkup
            <*> f ggg
            <*> f sl
            <*> f ntr
            <*> f c3
            <*> (denseV.vectorTraverse) f ht

deriving instance
  ( Show c, Show (LimitType c), Show e, Show (vc c), Show (ve e)
  ) => Show (Hairpin vc ve c e)

-- | Stacking helix contributions. Contains parameters for canonical stacks,
-- bulges, and interior loops.

data Stack ve c e = Stack
  { _stacking :: !(Dense ve (Z:.c:.c:.c:.c) e)
  }

stackE :: (VG.Vector ve a, VG.Vector we b) => Traversal (Stack ve c a) (Stack we c b) a b
{-# Inlinable stackE #-}
stackE f (Stack s) = Stack <$> (denseV.vectorTraverse) f s



-- * The full model.

data Turner2004 vc ve c e = Turner2004
  { _hairpin :: !(Hairpin vc ve c e)
    -- ^ Hairpin contributions
  , _stack   :: !(Stack ve c e)
  }
makeLenses ''Turner2004

-- | Traverse all energy fields with a polymorphic and energy-vector changing traversal.

turner2004E :: (VG.Vector ve a, VG.Vector we b) => Traversal (Turner2004 vc ve c a) (Turner2004 vc we c b) a b
{-# Inlinable turner2004E #-}
turner2004E f (Turner2004 hp s) = Turner2004 <$> hairpinE f hp <*> stackE f s



-- * Testing

test = Hairpin
  { _hairpinLength = VU.singleton (1::Int)
  , _hairpinMM = fromAssocs (ZZ:..aa:..aa:..aa:..aa) 2 []
  , _hairpinLookup = M.singleton (VU.singleton U) 3
  , _hairpinGGG = 4
  , _hairpinCslope = 5
  , _hairpinCintercept = 6
  , _hairpinC3 = 7
  , _hairpinTerm = fromAssocs (ZZ:..aa:..aa) 8 []
  }
  where aa= LtLetter A

