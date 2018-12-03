
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
--import Test.QuickCheck.Modifiers
--import Test.QuickCheck.Property
import Test.Tasty
import Test.Tasty.HUnit
--import Test.Tasty.QuickCheck
import Test.Tasty.TH
import Data.Maybe

import Data.PrimitiveArray ( (!), Z(..), (:.)(..) )
import Biobase.Primary.Nuc.RNA
import Biobase.Secondary.Vienna

import Biobase.Turner
import Biobase.Turner.Model.Vienna



---- * Tests of the ViennaRNA importer
--
---- | Make sure that we can parse successfully.
--
--case_Vienna_RNA_import = do
--  t ← viennaFromFile "./data/rna_turner2004.par"
--  assertBool "Import.Vienna parses successfully" $ isJust t
--
--
--
---- * Make sure that the tables are filled according to the layout of the
---- @.par@ file.
--
---- | The stack table.
--
--case_Vienna_RNA_table_stack = do
--  Just (e,p) ← viennaFromFile "./data/rna_turner2004.par"
--  assertEqual "G-C N_S" (-150) $ (e ^. stack) ! (Z:.GC:.NS)
--  assertEqual "G-C G-C" (-340) $ (e ^. stack) ! (Z:.GC:.GC)
--  assertEqual "G-C G-U" (-250) $ (e ^. stack) ! (Z:.GC:.GU)
--  assertEqual "A-U A-U" (-110) $ (e ^. stack) ! (Z:.AU:.AU)
--
--
--
---- * Make sure that the energy-calculating functions do the right thing.
---- They accept input characters in one order and do the correct lookup.
--
--case_Model_Vienna_eStack = do
--  Just (e,p) ← viennaFromFile "./data/rna_turner2004.par"
--  assertEqual "C A U G" (-210) $ eStack e C A U G
--  -- wc-sc-example.html
--  -- ~AGCGCU~AGCGCU~
--  -- ~[[[[[(~)]]]]]~
--  assertEqual "A G C U" (-210) $ eStack e A G C U
--  assertEqual "G C G C" (-340) $ eStack e G C G C
--  assertEqual "C G C G" (-240) $ eStack e C G C G
--  assertEqual "G C G C" (-340) $ eStack e G C G C
--  assertEqual "C U A G" (-210) $ eStack e C U A G

-- * Tests of the Turner importer

--case_Turner_RNA_Energy_import_ = do
--  t <- fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
--  return ()
--
--case_Turner_RNA_Energy_stack__ = do
--  t <- fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
--  assertEqual "G-G G-G stack" def    $ (t ^. stack) ! (Z:.G:.G:.G:.G)
--  assertEqual "G-C G-C stack" (-3.4) $ (t ^. stack) ! (Z:.G:.C:.G:.C)
--  assertEqual "G-C G-U stack" (-2.5) $ (t ^. stack) ! (Z:.G:.C:.G:.U)
----  assertEqual "A-U A-U stack" (-0.9) $ (t ^. stack) ! (Z:.A:.U:.A:.U)
--
--case_Turner_RNA_Energy_dangle3 = do
--  t <- fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
--  return ()
--
--case_Turner_RNA_Energy_dangle5 = do
--  t <- fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
--  return ()



main :: IO ()
main = $(defaultMainGenerator)

