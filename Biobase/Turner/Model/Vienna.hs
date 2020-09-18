
-- | This module provides a small number of functions that together implement the "Vienna-variant"
-- of the Turner energy model. This variant is being used in the ViennaRNA package and replicated
-- here for ease-of-use.
--
-- TODO benchmark variant separations, further variants where necessary.

module Biobase.Turner.Model.Vienna where

import           Data.Vector.Unboxed (Unbox)
import qualified Data.Map.Strict as MS
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Unboxed as VU

import           Data.PrimitiveArray
import           Algebra.Structure.Semiring

import           Biobase.Primary.Letter
import           Biobase.Primary.Nuc.RNA
import           Biobase.Secondary.Vienna
import           Biobase.Turner.Types
import           Biobase.Types.Energy



-- ** Stack variants: canonical stack, left and right bulge, interior loop.
-- Common function wrapping all four variants.

-- | Calculate the energy of a stack. The input order is the order of
-- characters on the tape. Say, @... C A ... U G ...@ leads to the call
-- @eStack C A U G@. The two outer characters @C - G@ are forming the new
-- pair on the stack, while the two inner characters @A - U@ form the old
-- pair on which to stack.
--
-- TODO check if @f@ is completely inlined. It is typically @\a b -> pairTable ! (Z:.a:.b)@
--
-- TODO @Letter u -> c@, since we can even be agnostic over what we index with ...

eStack
  :: ( Index c, VG.Vector ve e )
  => Stack ve c e
  -> c -> c -> c -> c
  -> e
{-# Inline eStack #-}
eStack Stack{..} l lp rp r =
  _stacking ! (Z:.l:.lp:.rp:.r)



-- | Generic helix, which calls the correct specialized function. This function calculates only the
-- helix energy of two pairs, without the continuation from a previous "weak" structure.

eHelix
  :: ( Index c, VG.Vector vc c, VG.Vector ve e )
  => Stack ve c e
  -> vc c -> vc c
  -> e
{-# Inline eHelix #-}
eHelix stack@Stack{..} ls rs
  | VG.length ls == 2, VG.length rs == 2
    = eStack stack (VG.unsafeHead ls) (VG.unsafeLast ls) (VG.unsafeHead rs) (VG.unsafeLast rs)


-- ** Hairpin

-- | The energy of a hairpin.
--
-- TODO check if @lenE@ is correct for hairpin >= 30.
--
-- TODO this is currently a somewhat simplified model not following the
-- @NNDB@ exactly.

eHairpin
  :: ( Semiring e, VG.Vector vc c, VG.Vector ve e, Index c, Ord c, Ord (vc c) )
  => Hairpin vc ve c e
  -> vc c
  -> e
{-# Inline eHairpin #-}
eHairpin Hairpin{..} xs
  -- disallow small loops
  | VG.length xs < 5 = zero
  -- look up special loops
  | Just e <- MS.lookup xs _hairpinLookup = e
  -- special case of loops with length 3
  | VG.length xs == 5 = mmE ⊗ lenE ⊗ term
  -- standard loops of no special case
  -- TODO missing special closures, as defined in the @NNDB@
  | otherwise = mmE ⊗ lenE ⊗ term
  where n = VG.length xs
        p = VG.unsafeIndex xs 0
        l = VG.unsafeIndex xs 1
        r = VG.unsafeIndex xs (n-1)
        q = VG.unsafeIndex xs n
        term = _hairpinTerm ! (Z:.p:.q)
        mmE  = _hairpinMM ! (Z:.p:.q:.l:.r)
        lenE = maybe lrgE
               id
             $ _hairpinLength VG.!? VG.length xs
        lrgE = error "fixme: lrgE"
--        lrgE = (_hairpinL VG.! 30) + (DekaG . round
--             . (*) (fromIntegral . getDekaG $ _largeLoop)
--             . log . (subtract 2) . fromIntegral $ VG.length xs)



-- ** Closure of a multibranched loop

--eMulti
--  :: ( Semiring e, VU.Unbox c )
--  => Multi c e
--  -> VU.Vector c
--  -> VU.Vector c
--  -> e
--{-# Inline eMulti #-}
--eMulti Multi{} ls rs = error "write me: eMulti"

