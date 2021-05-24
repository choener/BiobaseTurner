
-- | This module provides a small number of functions that together implement the "Vienna-variant"
-- of the Turner energy model. This variant is being used in the ViennaRNA package and replicated
-- here for ease-of-use.
--
-- TODO benchmark variant separations, further variants where necessary.

module Biobase.Turner.Model.Vienna where

import Data.Maybe (fromMaybe, isNothing)
import Data.Vector.Unboxed (Unbox)
import Debug.Trace (trace,traceShow)
import qualified Data.Map.Strict as MS
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Unboxed as VU
import Text.Printf (printf)

import Algebra.Structure.Semiring
import Data.PrimitiveArray as PA

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
-- Note: outer pair, inner pair, left and right unpaired

eIntLoop1x1
  :: (Index c, VG.Vector ve e)
  => IntLoop ve c e
  -> c -> c -> c -> c -> c -> c
  -> e
{-# Inline eIntLoop1x1 #-}
eIntLoop1x1 IntLoop{..} ol or il ir ul ur = -- outer, inner, unpaired, each left and right
  _intLoop1x1 ! (Z:.ol:.or:.il:.ir:.ul:.ur)

-- -- |
-- --
-- -- TODO fix up the order of characters!
-- 
eIntLoop1x2, eIntLoop2x1
  :: (Index c, VG.Vector ve e)
  => IntLoop ve c e
  -> c -> c -> c -> c -> c -> c -> c
  -> e
{-# Inline eIntLoop1x2 #-}
{-# Inline eIntLoop2x1 #-}
eIntLoop1x2 IntLoop{..} ol or il ir ul m2r m1r =
  _intLoop1x2 ! (Z:.ol:.or:.il:.ir:.ul:.m2r:.m1r)
eIntLoop2x1 IntLoop{..} ol or il ir u1p u2p ur =
  _intLoop1x2 ! (Z:.ol:.or:.il:.ir:.ur:.u1p:.u2p)

-- |
--
-- TODO fix up the order of characters!

eIntLoop2x2
  :: (Index c, VG.Vector ve e)
  => IntLoop ve c e
  -> c -> c -> c -> c -> c -> c -> c -> c
  -> e
{-# Inline eIntLoop2x2 #-}
eIntLoop2x2 IntLoop{..} lsO rsO rsI lsI lp1 lp2 rm2 rm1 =
  _intLoop2x2 ! (Z:.lsO:.rsO:.rsI:.lsI:.lp1:.lp2:.rm2:.rm1)

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
  :: ( Index c, VG.Vector vc c, VG.Vector ve e, Semiring e, Show (vc c), Ord e, Show c, Show e, Read c, Eq c )
  => Stack ve c e
  -> IntLoop ve c e
  -> vc c -> c
  -> c -> vc c
  -> e
{-# Inline scoreInteriorLoop #-}
scoreInteriorLoop stack@Stack{..} intloop@IntLoop{..} ls lsI rsI rs -- ls, lsInner, rsInner, rs
  -- Canonical stack 2x2
  | lsLen==1, rsLen==1 = eStack stack lsO rsO rsI lsI -- note outer, inner pair, both 5'--3'
  -- Stack with 1nt slippage: NOTE head/last switched, because of the "free" nt. In principle, we
  -- can unify canonical and slippage stacks
  | short==1&&long==2 = eStack stack lsO rsO rsI lsI ⊗ _bulgeL!(Z:.1)
  -- 1x1 interior loop
  | short==2, long==2 = let e = eIntLoop1x1 intloop lsO rsO rsI lsI lp1 rm1
                        in  e
  -- 2x1 interior loops
  | lsLen==2, rsLen==3 = eIntLoop1x2 intloop lsO rsO rsI lsI lp1 rm2 rm1
  | lsLen==3, rsLen==2 = eIntLoop2x1 intloop lsO rsO rsI lsI lp1 lp2 rm1
  -- 2x2 interior loop
  | lsLen==3, rsLen==3 = eIntLoop2x2 intloop lsO rsO rsI lsI lp1 lp2 rm2 rm1
--  -- 2x3 interior loops
--  | min lls rrs==3, max lls rrs==4 = _intLoop2x3!(Z:.lo0:.ro0:.lo1:.ro1)
--                                   ⊗ _intLoop2x3!(Z:.ri0:.li0:.ri1:.li1)
--                                   ⊗ _intLoopL!(Z:.5) ⊗ _intLoopNinio
  -- bulges to the right or left
--  |  (lls==1 && rrs>3 && rrs<31)
--  || (lls>3 && lls<31 && rrs==1) = _bulgeAU!(Z:.lo0:.ro0)
--                                 ⊗ _bulgeAU!(Z:.ri0:.li0)
--                                 ⊗ _bulgeL!(Z:.max (lls-2) (rrs-2))
  -- 1xn loop to the right
  |  (short==2 && long>3 && long<31) = let o = _intLoop1xn!(Z:.lsO:.rsO:.lp1:.rm1)
                                           i = _intLoop1xn!(Z:.rsI:.lsI:.rp1:.lm1)
                                           l = _intLoopL!(Z:.long + short - 2)
                                           n = one -- min _intLoopMaxNinio ((rrs-2) `nTimes` _intLoopNinio)
                                           e = o ⊗ i ⊗ l -- ⊗ n
                                       in  e
                                      -- if (lsO == read "G" && rsO == read "C" && lsI == read "C" && rsI == read "G")
                                      --     then traceShow (lsO,rsO,lp1,rm1,rsI,lsI,rp1,lm1,o,i,l,e) e
                                      --     else e
--  -- generic interior loops
  | short>3, lsLen+rsLen<34 = let o = _intLoopMM!(Z:.lsO:.rsO:.lp1:.rm1)
                                  i = _intLoopMM!(Z:.rsI:.lsI:.rp1:.lm1)
                                  l = _intLoopL!(Z:.long + short - 2)
                                  n = one -- min _intLoopMaxNinio (abs (lls-rrs) `nTimes` _intLoopNinio)
                                  e = o ⊗ i ⊗ l ⊗ n
                              in  e
--  -- for every other loop, we do not have energy parameters
--  | otherwise = trace (printf "scoreInteriorLoop too large: ls:%s rs:%s\n" (show ls) (show rs)) zero
  | otherwise = zero
  where
    lsLen = VG.length ls
    rsLen = VG.length rs
    short = min lsLen rsLen
    long  = max lsLen rsLen
    -- the outer pair. As @lsO ... inner part ... rsO@
    lsO = VG.unsafeHead ls
    rsO = VG.unsafeLast rs
    lp1 = VG.unsafeIndex ls 1
    lp2 = VG.unsafeIndex ls 2
    rm1 = VG.unsafeIndex rs (rsLen-2) -- second to last!
    rm2 = VG.unsafeIndex rs (rsLen-3)
    --
    lm1 = VG.unsafeIndex ls (lsLen-1) -- note that this is the last unpaired on the ls side
    rp1 = VG.unsafeIndex rs 0 -- note that this is the "last" unpaired on the rs side

--    -- define nucleotides for the left @ls@ and right @rs@ input string. @o@ defines the outside
--    -- character, leftmost in @l@, rightmost in @r@. The index is the offset from this character.
--    -- Hence @lo0@ is the leftmost character of the left string, while @ro0@ the rightmost character
--    -- of the right string.
--    lo0 = ls VG.! 0;       li0 = lsI
--    lo1 = ls VG.! 1;       li1 = ls VG.! (lls-1)
--    lo2 = ls VG.! 2;       li2 = ls VG.! (lls-2)
--    ri0 = rs VG.! (rrs-1); ro0 = rs VG.! 0
--    ri1 = rs VG.! (rrs-2); ro1 = rs VG.! 1
--    ri2 = rs VG.! (rrs-3); ro2 = rs VG.! 2
--    -- do not evaluate strictly, some branches are possibly illegal!
--    ls1 = VG.unsafeIndex ls 1; ls2 = VG.unsafeIndex ls 2
--    rs1 = VG.unsafeIndex rs 1; rs2 = VG.unsafeIndex rs 2


-- ** The hairpin loop

-- | The energy of a hairpin.
--
-- HP(xs) = mismatch(x_0, x_n, x_1, x_{n-1}) + length(n-2)
--
-- TODO incomplete

scoreHairpin
  :: ( Semiring e, VG.Vector vc c, VG.Vector ve e, Index c, Ord c, Ord (vc c), Show e, Show c, Show (vc c) )
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
  | VG.length xs == 5 = lenE ⊗ wobble
  -- standard loops of no special case
  -- TODO missing special closures, as defined in the @NNDB@
  | otherwise = _hairpinMM!(Z:.lo0:.ro0:.lo1:.ro1) ⊗ lenE
  where n = VG.length xs
        lo0 = VG.unsafeIndex xs 0;     lo1 = VG.unsafeIndex xs 1
        ro0 = VG.unsafeIndex xs (n-1); ro1 = VG.unsafeIndex xs (n-2)
        lenE = fromMaybe lrgE $ _hairpinLength !? (Z:.n-2)
        llpE = (floor $ log (fromIntegral n-2) / 30) `nTimes` _largeLoop
        lrgE = (_hairpinLength!(Z:.30)) ⊗ llpE
        wobble = _hairpinPenalty ! (Z:.lo0:.ro0)

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

-- ** Score external loops

-- | If the external loop has characters to left or right, score with a mismatch system.

scoreExteriorLoop
  :: ( Semiring e, VU.Unbox c, Index c, VG.Vector ve e, Show e, Show c )
  => Exterior ve c e
  -> Maybe c -> Maybe c
  -> Maybe c -> Maybe c
  -> e
{-# Inline scoreExteriorLoop #-}
scoreExteriorLoop Exterior{..} mlo mli mri mro
  -- mismatch exterior
  -- TODO the end penalty follows VRNA, but seems dubious
  | Just lo <- mlo, Just li <- mli
  , Just ri <- mri, Just ro <- mro = let e = _mismatchExterior ! (Z:.li:.ri:.lo:.ro)
                                           ⊗ _endPenalty ! (Z:.li:.ri)
                                     in  e
  -- dangle5
  -- TODO the end penalty follows VRNA, but seems dubious
  | Just lo <- mlo, Nothing <- mro
  , Just li <- mli, Just ri <- mri = let e = _dangle5 ! (Z:.li:.ri:.lo)
                                           ⊗ _endPenalty ! (Z:.li:.ri)
                                     in  e
  -- dangle3
  -- TODO the end penalty follows VRNA, but seems dubious
  | Nothing <- mlo, Just ro <- mro
  , Just li <- mli, Just ri <- mri = let e = _dangle3 ! (Z:.li:.ri:.ro)
                                           ⊗ _endPenalty ! (Z:.li:.ri)
                                     in  e
  | Nothing <- mlo, Nothing <- mro
  , Just li <- mli, Just ri <- mri = let e = _endPenalty ! (Z:.li:.ri)
                                     in  e
  | otherwise = one

