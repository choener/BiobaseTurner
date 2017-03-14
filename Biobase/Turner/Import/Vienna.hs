
-- | Import a @ViennaRNA 2004@ energy parameter file.

module Biobase.Turner.Import.Vienna where

import Text.Trifecta as TT
import Data.Text (Text,pack)
import Data.PrimitiveArray as PA hiding (map)
import Control.Monad (guard)
import Control.Applicative ( (<|>) )

import Biobase.Turner.Types
import Biobase.Turner.Import.Turner (minPP,maxPP,minPBB,maxPBB)
import Biobase.Types.Energy



fromFile :: FilePath -> IO (Maybe Vienna2004)
fromFile = TT.parseFromFile pVienna

pVienna :: Parser Vienna2004
pVienna = v2Header *> spaces *> (blocksToTurner <$> some (block <* whiteSpace)) <* v2End <* eof
  where v2Header = text "## RNAfold parameter file v2.0" <?> ".par Header"
        v2End = text "#END"

blocksToTurner = error "blocksToTurner"

-- | Blocks have different arrays. The block parser will switch into the
-- corresponding parser where necessary.

data Block
  = BlockPP Text EE (Arr PP)
  | BlockPBB Text EE (Arr PNN)

data EE = Entropy | Enthalpie

type Arr i = Unboxed i DeltaDekaGibbs

block :: Parser Block
block = blockPP <|> blockPBB

blockPP :: Parser Block
blockPP = (\(t,e) vs -> BlockPP t e vs) <$> (try $ headerParser blockHeader) <*> (ppArray >>= convertPP) <?> "blockPP"
  where convertPP :: [Int] -> Parser (Arr PP)
        convertPP vs = do
          guard (length vs == 7*7) <?> "blockPP with wrong number of parsed integers (" ++ show (length vs) ++ ")"
          return $ fromAssocs minPP maxPP (DekaG 999999) [error "blockPP"]
        blockHeader = ["stack"]

blockPBB :: Parser Block
blockPBB = (\(t,e) vs -> BlockPBB t e vs) <$> (try $ headerParser blockHeader) <*> (pnnArray >>= convertPBB) <?> "blockPBB"
  where convertPBB vs = do
          guard (length vs == 7*5*5) <?> "blockPBB with wrong number of parsed integers (" ++ show (length vs) ++ ")"
          return $ fromAssocs minPBB maxPBB (DekaG 999999) [error "blockPBB"]
        blockHeader = ["mismatch_hairpin", "mismatch_interior", "mismatch_interior_1n", "mismatch_interior_23", "mismatch_multi", "mismatch_exterior"]

-- | Try different headers, succeed only if the full line can be parsed.
-- Required because prefixes of headers overlap.

headerParser :: [Text] -> Parser (Text,EE)
headerParser ps = choice [try $ (,) <$ text "#" <* whiteSpace <*> text p <*> pEE | p <- ps]
  where pEE = (Entropy <$ newline) <|> (Enthalpie <$ text "_enthalpies" <* newline)

ppArray :: Parser [Int]
ppArray = do
  headerComment <- commentLine
  valueLines <- count 7 valueLine
  return $ concat valueLines

commentLine = optional $ between (text "/*") (text "*/") (spaces *> some (letter <|> oneOf ",") `sepEndBy` spaces) <?> "commentLine"

valueLine :: Parser [Int]
valueLine = spaces *> (fromInteger <$> integer) `sepEndBy` spaces <* commentLine <* newline

pnnArray :: Parser [Int]
pnnArray = do
  headerComment <- commentLine
  valueLines <- count (7*5) valueLine
  return $ concat valueLines


test = fromFile "data/rna_turner2004.par"

