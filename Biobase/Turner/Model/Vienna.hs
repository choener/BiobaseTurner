
-- | This module provides a small number of functions that together implement the "Vienna-variant"
-- of the Turner energy model. This variant is being used in the ViennaRNA package and replicated
-- here for ease-of-use.
--
-- TODO benchmark variant separations, further variants where necessary.

module Biobase.Turner.Model.Vienna where

import Data.Maybe (fromMaybe)
import Data.Vector.Unboxed (Unbox)
import Debug.Trace (trace,traceShow)
import qualified Data.Map.Strict as MS
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Unboxed as VU
import Text.Printf (printf)

import Algebra.Structure.Semiring
import Data.PrimitiveArray

import Biobase.Primary.Letter
import Biobase.Primary.Nuc.RNA
import Biobase.Secondary.Vienna
import Biobase.Turner.Types
import Biobase.Types.Energy



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
  :: ( Index c, VG.Vector ve e, Show c, Show e )
  => Stack ve c e
  -> c -> c -> c -> c
  -> e
{-# Inline eStack #-}
eStack Stack{..} l lp rp r =
  let e = _stacking ! (Z:.l:.lp:.rp:.r)
  in  e -- traceShow (l,lp,rp,r,e) e

-- |
--
-- TODO fix up the order of characters!

eIntLoop1x1
  :: (Index c, VG.Vector ve e)
  => IntLoop ve c e
  -> c -> c -> c -> c -> c -> c
  -> e
{-# Inline eIntLoop1x1 #-}
eIntLoop1x1 IntLoop{..} l lm lp rp rm r =
  _intLoop1x1 ! (Z:.l:.lm:.lp:.rp:.rm:.r)

-- -- |
-- --
-- -- TODO fix up the order of characters!
-- 
-- eIntLoop1x2, eIntLoop2x1
--   :: (Index c, VG.Vector ve e)
--   => IntLoop ve c e
--   -> c -> c -> c -> c -> c -> c -> c
--   -> e
-- {-# Inline eIntLoop1x2 #-}
-- {-# Inline eIntLoop2x1 #-}
-- eIntLoop1x2 IntLoop{..} l l1 lp rp r1 r2 r =
--   _intLoop1x2 ! (Z:.l:.l1:.lp:.rp:.r1:.r2:.r)
-- eIntLoop2x1 IntLoop{..} l l1 l2 lp rp r1 r =
--   _intLoop1x2 ! (Z:.l:.l1:.l2:.lp:.rp:.r1:.r)
-- 
-- -- |
-- --
-- -- TODO fix up the order of characters!
-- 
-- eIntLoop2x2
--   :: (Index c, VG.Vector ve e)
--   => IntLoop ve c e
--   -> c -> c -> c -> c -> c -> c -> c -> c
--   -> e
-- {-# Inline eIntLoop2x2 #-}
-- eIntLoop2x2 IntLoop{..} l l1 l2 lp rp r1 r2 r =
--   _intLoop2x2 ! (Z:.l:.l1:.l2:.lp:.rp:.r1:.r2:.r)

-- | Generic interior loop, which calls the correct specialized function. This function calculates
-- only the loop energy of two pairs, without the continuation from a previous "weak" structure. We
-- have
-- @
-- l:ls ... rs:r
-- [l0,l1,...,r1,r0]
--  (  ...()...  )
-- @
--
-- 2x3 loops
-- @
-- 1234567890
-- (..(-)...)
-- @
--
-- TODO Expand to a model of "unlimited" size, which will allow handing in inputs of any size.

scoreInteriorLoop
  :: ( Index c, VG.Vector vc c, VG.Vector ve e, Semiring e, Show (vc c), Ord e, Show c, Show e )
  => Stack ve c e
  -> IntLoop ve c e
  -> vc c -> c
  -> c -> vc c
  -> e
{-# Inline scoreInteriorLoop #-}
scoreInteriorLoop stack@Stack{..} intloop@IntLoop{..} ls lsR rsL rs
  -- Canonical stack 2x2
  | lls==1, rrs==1 = eStack stack (VG.unsafeLast ls) (VG.unsafeHead rs) rsL lsR -- note outer, inner pair, both 5'--3'
  -- Stack with 1nt slippage: NOTE head/last switched, because of the "free" nt. In principle, we
  -- can unify canonical and slippage stacks
  | lls==1&&rrs==2 || lls==2&&rrs==1 = eStack stack (VG.unsafeHead ls) (VG.unsafeLast rs) rsL lsR -- ⊗ _bulgeL!(Z:.1)
--  -- 1x1 interior loop
--  | lls==2, rrs==2 = eIntLoop1x1 intloop lo0 ro0 ri0 li0 lo1 ro1
--  -- 2x1 interior loops
--  | lls==2, rrs==3 = eIntLoop1x2 intloop lo0 ro0 ri0 li0 lo1 ri1 ri2
--  | lls==3, rrs==2 = eIntLoop2x1 intloop lo0 ro0 ri0 li0 ro1 li2 li1
--  -- 2x2 interior loop
--  | lls==3, rrs==3 = eIntLoop2x2 intloop lo0 ro0 ri0 li0 lo1 lo2 ri2 ri1
--  -- 2x3 interior loops
--  | min lls rrs==3, max lls rrs==4 = _intLoop2x3!(Z:.lo0:.ro0:.lo1:.ro1)
--                                   ⊗ _intLoop2x3!(Z:.ri0:.li0:.ri1:.li1)
--                                   ⊗ _intLoopL!(Z:.5) ⊗ _intLoopNinio
--  -- bulge to the right
--  |  (lls==1 && rrs>3 && rrs<31)
--  || (lls>3 && lls<31 && rrs==1) = _bulgeAU!(Z:.lo0:.ro0)
--                                 ⊗ _bulgeAU!(Z:.ri0:.li0)
--                                 ⊗ _bulgeL!(Z:.max (lls-2) (rrs-2))
--  -- 1xn loop to the right
--  |  (lls==2 && rrs>3 && rrs<31)
--  || (rrs==2 && lls>3 && lls<31) = _intLoop1xn!(Z:.lo0:.ro0:.lo1:.ro1)
--                                 ⊗ _intLoop1xn!(Z:.ri0:.li0:.ri1:.li1)
--                                 ⊗ _intLoopL!(Z:.rrs-2+1)
--                                 ⊗ min _intLoopMaxNinio ((rrs-2) `nTimes` _intLoopNinio)
--  -- generic interior loops
--  | lls>3, rrs>3, lls+rrs>8, lls+rrs<34 = _intLoopMM!(Z:.lo0:.ro0:.lo1:.ro1)
--                                        ⊗ _intLoopMM!(Z:.ri0:.li0:.ri1:.li1)
--                                        ⊗ _intLoopL!(Z:.lls+rrs-4)
--                                        ⊗ min _intLoopMaxNinio (abs (lls-rrs) `nTimes` _intLoopNinio)
--  -- for every other loop, we do not have energy parameters
--  | otherwise = trace (printf "scoreInteriorLoop too large: ls:%s rs:%s\n" (show ls) (show rs)) zero
  | otherwise = zero
  where
    lls = VG.length ls
    rrs = VG.length rs
    -- define nucleotides for the left @ls@ and right @rs@ input string. @o@ defines the outside
    -- character, leftmost in @l@, rightmost in @r@. The index is the offset from this character.
    -- Hence @lo0@ is the leftmost character of the left string, while @ro0@ the rightmost character
    -- of the right string.
    lo0 = ls VG.! 0;       li0 = lsR
    lo1 = ls VG.! 1;       li1 = ls VG.! (lls-1)
    lo2 = ls VG.! 2;       li2 = ls VG.! (lls-2)
    ri0 = rs VG.! (rrs-1); ro0 = rs VG.! 0
    ri1 = rs VG.! (rrs-2); ro1 = rs VG.! 1
    ri2 = rs VG.! (rrs-3); ro2 = rs VG.! 2
    -- do not evaluate strictly, some branches are possibly illegal!
    ls1 = VG.unsafeIndex ls 1; ls2 = VG.unsafeIndex ls 2
    rs1 = VG.unsafeIndex rs 1; rs2 = VG.unsafeIndex rs 2


-- ** The hairpin loop

-- | The energy of a hairpin.
--
-- HP(xs) = mismatch(x_0, x_n, x_1, x_{n-1}) + length(n-2)
--
-- TODO incomplete

scoreHairpin
  :: ( Semiring e, VG.Vector vc c, VG.Vector ve e, Index c, Ord c, Ord (vc c) )
  => Hairpin vc ve c e
  -> vc c
  -> e
{-# Inline scoreHairpin #-}
scoreHairpin Hairpin{..} xs
  -- look up special loops
  | Just e <- MS.lookup xs _hairpinLookup = e
  -- disallow small loops
  | VG.length xs < 5 = zero
  -- special case of loops with length 3
  | VG.length xs == 5 = lenE -- ⊗ _hairpinAU ! (Z:.lo0:.ro0)
  -- standard loops of no special case
  -- TODO missing special closures, as defined in the @NNDB@
  | otherwise = _hairpinMM!(Z:.lo0:.ro0:.lo1:.ro1) ⊗ lenE
  where n = VG.length xs
        lo0 = VG.unsafeIndex xs 0;     lo1 = VG.unsafeIndex xs 1
        ro0 = VG.unsafeIndex xs (n-1); ro1 = VG.unsafeIndex xs (n-2)
        lenE = fromMaybe lrgE $ _hairpinLength !? (Z:.n-2)
        llpE = (floor $ log (fromIntegral n-2) / 30) `nTimes` _largeLoop
        lrgE = (_hairpinLength!(Z:.30)) ⊗ llpE

-- ** Closure of a multibranched loop

scoreMultiLoop
  :: ( Semiring e, VU.Unbox c, Index c, VG.Vector ve e )
  => Multi ve c e
  -> c -> c
  -> c -> c
  -> e
{-# Inline scoreMultiLoop #-}
scoreMultiLoop Multi{..} l lp rp r =
  let e = _mismatchMulti ! (Z:.l:.lp:.rp:.r)
  in  e -- traceShow (l,lp,rp,r,e) e

