
-- | Import a @ViennaRNA 2004@ energy parameter file.

module Biobase.Turner.Import.Vienna where

import           Control.Applicative ( (<|>) )
import           Control.Monad (guard)
import           Data.PrimitiveArray as PA hiding (map,fromList)
import           Data.Text (Text,pack)
import           Data.Vector.Unboxed (Vector,fromList)
import           Text.Trifecta as TT
import           Data.Map.Strict (Map)
import qualified Data.Map.Strict as M

import Biobase.Primary (Primary,RNA,primary)
import Biobase.Turner.Import.Turner (minPP,maxPP,minPBB,maxPBB,minPB,maxPB,minPPBB,maxPPBB,minPPBBB,maxPPBBB,minPPBBBB,maxPPBBBB)
import Biobase.Turner.Types
import Biobase.Types.Energy



fromFile :: FilePath -> IO (Maybe Vienna2004)
fromFile = TT.parseFromFile pVienna

pVienna :: Parser Vienna2004
pVienna = v2Header *> spaces *> (blocksToTurner <$> some (block <* whiteSpace)) <* v2End <* eof
  where v2Header = text "## RNAfold parameter file v2.0" <?> ".par Header"
        v2End = text "#END" <* newline

blocksToTurner = error "blocksToTurner"

-- | Blocks have different arrays. The block parser will switch into the
-- corresponding parser where necessary.

data Block
  = BlockPP     Text EE (Arr PP)
  | BlockPBB    Text EE (Arr PNN)
  | BlockPB     Text EE (Arr PN)
  | BlockPPBB   Text EE (Arr PPNN)
  | BlockPPBBB  Text EE (Arr PPNNN)
  | BlockPPBBBB Text EE (Arr PPNNNN)
  | BlockLinear Text EE (Vector DeltaDekaGibbs)
  | BlockML     Text DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs
  | BlockNinio  Text DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs
  | BlockMisc   Text DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs
  | BlockLoops  Text    (Map (Primary RNA) (DeltaDekaGibbs,DeltaDekaGibbs))

data EE = Entropy | Enthalpie

type Arr i = Unboxed i DeltaDekaGibbs

block :: Parser Block
block =   blockPP <|> blockPBB <|> blockPB <|> blockPPBB <|> blockPPBBB <|> blockPPBBBB
      <|> blockLinear <|> blockML <|> blockNinio <|> blockMisc <|> blockLoops

blockPP :: Parser Block
blockPP = wrapBlock BlockPP blockHeader convertPP
  where convertPP vs = do
          let l = length vs
          guard (l == 7*7) <?> "expected 7*7=49 entries for blockPP but got " ++ show l
          return $ fromAssocs minPP maxPP (DekaG 999999) [error "blockPP"]
        blockHeader = ["stack"]

blockPBB :: Parser Block
blockPBB = wrapBlock BlockPBB blockHeader convertPBB
  where convertPBB vs = do
          let l = length vs
          guard (l == 7*5*5) <?> "expected 7*5*5=175 entries for blockPBB but got " ++ show l
          return $ fromAssocs minPBB maxPBB (DekaG 999999) [error "blockPBB"]
        blockHeader = ["mismatch_hairpin", "mismatch_interior", "mismatch_interior_1n", "mismatch_interior_23", "mismatch_multi", "mismatch_exterior"]

blockPB :: Parser Block
blockPB = wrapBlock BlockPB blockHeader convertPB
  where convertPB vs = do
          let l = length vs
          guard (l == 7*5) <?> "7*5=35 entries for blockPB but got " ++ show l
          return $ fromAssocs minPB maxPB (DekaG 999999) [error "blockPB"]
        blockHeader = ["dangle5", "dangle3"]

blockPPBB :: Parser Block
blockPPBB = wrapBlock BlockPPBB blockHeader convertPPBB
  where convertPPBB vs = do
          let l = length vs
          guard (l == 7*7*5*5) <?> "7*7*5*5=1225 entries for blockPPBB but got " ++ show l
          return $ fromAssocs minPPBB maxPPBB (DekaG 999999) [error "blockPPBB"]
        blockHeader = ["int11"]

blockPPBBB :: Parser Block
blockPPBBB = wrapBlock BlockPPBBB blockHeader convertPPBBB
  where convertPPBBB vs = do
          let l = length vs
          guard (l == 7*7*5*5*5) <?> "7*7*5*5*5=6125 entries for blockPPBBB but got " ++ show l
          return $ fromAssocs minPPBBB maxPPBBB (DekaG 999999) [error "blockPPBBB"]
        blockHeader = ["int21"]

-- | Parses the big @int22@ lists. Excludes special characters and only
-- includes standard nucleotides / pairs.

blockPPBBBB :: Parser Block
blockPPBBBB = wrapBlock BlockPPBBBB blockHeader convertPPBBBB
  where convertPPBBBB vs = do
          let l = length vs
          guard (l == 6*6*4*4*4*4) <?> "6*6*4*4*4*4=30625 entries for blockPPBBBB but got " ++ show l
          return $ fromAssocs minPPBBBB maxPPBBBB (DekaG 999999) [error "blockPPBBBB"]
        blockHeader = ["int22"]

blockLinear :: Parser Block
blockLinear = wrapBlock BlockLinear blockHeader convertLinear
  where convertLinear vs = do
          let l = length vs
          guard (l == 31) <?> "31 entries for blockLinear but got " ++ show l
          return $ fromList vs
        blockHeader = ["hairpin", "bulge", "interior"]

blockML :: Parser Block
blockML = BlockML <$> (try $ id <$ text "# " <*> text "ML_params" <* newline)
                  <*> ddg <*> ddg <*> ddg <*> ddg <*> ddg <*> ddg

blockNinio :: Parser Block
blockNinio = BlockNinio <$> (try $ id <$ text "# " <*> text "NINIO" <* newline)
                        <*> ddg <*> ddg <*> ddg

blockMisc :: Parser Block
blockMisc = BlockMisc <$> (try $ id <$ text "# " <*> text "Misc" <* newline)
                      <*> ddg <*> ddg <*> ddg <*> ddg

blockLoops :: Parser Block
blockLoops = choice [BlockLoops <$> (try $ id <$ text "# " <*> text lp <* newline) <*> loopLines | lp <- loops]
  where loopLines = M.fromList <$> some triplet -- `sepBy` newline
        triplet = (\k v w -> (k,(v,w))) <$> (primary <$> some (oneOf "ACGU") <* spaces)
                                        <*> ddg <*> ddg
        loops = ["Hexaloops", "Tetraloops", "Triloops"]

ddg :: Parser DeltaDekaGibbs
ddg = (DekaG . fromIntegral) <$ spaces <*> integer

wrapBlock :: (Text -> EE -> arr -> Block) -> [Text] -> ([DeltaDekaGibbs] -> Parser arr) -> Parser Block
wrapBlock block blockHeader convert
  =   (\(t,e) vs -> block t e vs) <$> headerParser blockHeader <*> (arrayLines >>= convert)
  <?> "wrapBlock"

-- | Try different headers, succeed only if the full line can be parsed.
-- Required because prefixes of headers overlap.

headerParser :: [Text] -> Parser (Text,EE)
headerParser ps = choice [try $ (,) <$ text "#" <* whiteSpace <*> text p <*> pEE | p <- ps]
  where pEE = (Entropy <$ newline) <|> (Enthalpie <$ text "_enthalpies" <* newline)

arrayLines :: Parser [DeltaDekaGibbs]
arrayLines = do
  headerComment <- commentLine
  valueLines <- valueLine `sepEndBy` newline
  return . map DekaG $ concat valueLines

commentLine = optional $ between (text "/*") (text "*/") (spaces *> some (letter <|> oneOf ",") `sepEndBy` spaces) <?> "commentLine"

valueLine :: Parser [Int]
valueLine = spaces *> (fromInteger <$> integer <|> 999999 <$ text "INF") `sepEndBy` spaces <* commentLine -- <* newline


test = fromFile "data/rna_turner2004.par"

