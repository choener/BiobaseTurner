
module Biobase.Turner.Model.Vienna where

import qualified Data.Vector.Generic as VG
import qualified Data.Map.Strict as MS

import           Data.PrimitiveArray

import           Biobase.Primary.Letter
import           Biobase.Primary.Nuc.RNA
import           Biobase.Secondary.Vienna
import           Biobase.Turner.Types
import           Biobase.Types.Energy



-- | Calculate the energy of a stack. The input order is the order of
-- characters on the tape. Say, @... C A ... U G ...@ leads to the call
-- @eStack C A U G@. The two outer characters @C - G@ are forming the new
-- pair on the stack, while the two inner characters @A - U@ form the old
-- pair on which to stack.

eStack ∷ Vienna2004 → Letter RNA → Letter RNA → Letter RNA → Letter RNA → DeltaDekaGibbs
eStack Turner2004Model{..} l lp rp r =
  let p = viennaPairTable ! (Z:.l:.r)
      q = viennaPairTable ! (Z:.rp:.lp)
  in  _stack ! (Z:.p:.q)
{-# Inline eStack #-}

-- | The energy of a hairpin.
--
-- TODO check if @lenE@ is correct for hairpin >= 30.
--
-- TODO this is currently a somewhat simplified model not following the
-- @NNDB@ exactly.

eHairpin ∷ Vienna2004 → Letter RNA → Vector (Letter RNA) → Letter RNA → DeltaDekaGibbs
eHairpin Turner2004Model{..} l us r
  -- disallow small loops
  | VG.length us < 3
  = DekaG 999999
  -- look up special loops
  | Just e ← MS.lookup lusr _hairpinLookup
  = e
  -- special case of loops with length 3
  | VG.length us == 3
  = mmE + lenE + tAUE
  -- standard loops of no special case
  -- TODO missing special closures, as defined in the @NNDB@
  | otherwise
  = mmE + lenE + tAUE
  where lusr = VG.snoc (VG.cons ll us) rr
        allC = VG.all (==C) us
        ll   = VG.unsafeHead us
        rr   = VG.unsafeLast us
        p    = viennaPairTable ! (Z:.l:.r)
        mmE  = _hairpinMM ! (Z:.p:.ll:.rr)
        lenE = maybe lrgE
               id
             $ _hairpinL VG.!? VG.length us
        lrgE = (_hairpinL VG.! 30) + (DekaG . round
             . (*) (fromIntegral . getDekaG $ _largeLoop)
             . log . fromIntegral $ VG.length us)
        tAUE = if p == AU then _termAU else 0
{-# Inline eHairpin #-}

