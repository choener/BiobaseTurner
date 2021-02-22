
-- | Default @Turner RNA 2004@ energy parameters parsed from
-- @data/rna_turner2004.par@.

module Biobase.Turner.Default.RNA2004
  ( module Biobase.Turner.Default.RNA2004
  , Vienna2004
  ) where

import Data.ByteString
import Data.FileEmbed

import Biobase.Turner.Types
import Biobase.Turner.Import.Vienna



-- | The raw parameter file embedded.

rawTurner2004 :: ByteString
rawTurner2004 = $(makeRelativeToProject "data/rna_turner2004.par" >>= embedFile)

-- | Parameters parsed into the entropy and enthalpy terms.

turner2004 :: (Vienna2004, Vienna2004)
turner2004 = case fromByteString rawTurner2004 of
  Left err -> error err
  Right xs -> xs

