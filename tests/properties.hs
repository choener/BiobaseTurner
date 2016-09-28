
-- | This property test shall read in Turner files and then check
-- individual contributions that are known to be right.
--
-- TODO Known energies for known structures will have to be tested in, for
-- example, the ViennaRNA-bindings; or our own Haskell package.
--
-- NOTE we run @fromDir@ quite often, but want to have each test separated.

module Main where

import Control.Lens
import Data.Default
import Test.QuickCheck.Modifiers
import Test.QuickCheck.Property
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck
import Test.Tasty.TH

import Data.PrimitiveArray ( (!), Z(..), (:.)(..) )
import Biobase.Primary.Nuc.RNA

import Biobase.Turner



case_Turner_RNA_Energy_import_ = do
  t <- fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
  return ()

case_Turner_RNA_Energy_stack__ = do
  t <- fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
  assertEqual "G-G G-G stack" def    $ (t ^. stack) ! (Z:.G:.G:.G:.G)
  assertEqual "G-C G-C stack" (-3.4) $ (t ^. stack) ! (Z:.G:.C:.G:.C)
  assertEqual "G-C G-U stack" (-2.5) $ (t ^. stack) ! (Z:.G:.C:.G:.U)
--  assertEqual "A-U A-U stack" (-0.9) $ (t ^. stack) ! (Z:.A:.U:.A:.U)

case_Turner_RNA_Energy_dangle3 = do
  t <- fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
  return ()

case_Turner_RNA_Energy_dangle5 = do
  t <- fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
  return ()



main :: IO ()
main = $(defaultMainGenerator)

