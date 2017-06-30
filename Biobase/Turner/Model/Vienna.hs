
module Biobase.Turner.Model.Vienna where

import Data.PrimitiveArray

import Biobase.Secondary.Vienna
import Biobase.Primary.Nuc.RNA
import Biobase.Primary.Letter
import Biobase.Types.Energy
import Biobase.Turner.Types



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

