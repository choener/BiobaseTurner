
-- | The 'Turner2004' data structure reflects the RNA (and DNA) energy
-- parameters known as the Turner 2004 data set.
--
-- In general, have a look here:
-- <http://rna.urmc.rochester.edu/NNDB/turner04/index.html> where parameters
-- are explained.

module Biobase.Turner
  ( module Biobase.Turner.Types
--  , turnerFromDir
--  , viennaFromFile
  ) where

import           Biobase.Turner.Types
import qualified Biobase.Turner.Import.Turner as T
import qualified Biobase.Turner.Import.Vienna as V

-- turnerFromDir ∷ String → String → String → IO Turner2004
-- turnerFromDir = T.fromDir
-- 
-- -- | Import a Vienna energy file. In case @Nothing@ is returned, the errors
-- -- are on the console and we should quit the program.
-- 
-- viennaFromFile ∷ String → IO (Maybe (Vienna2004, Vienna2004))
-- viennaFromFile = V.fromFile
-- 
