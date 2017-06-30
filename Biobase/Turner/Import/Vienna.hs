
-- | Import a @ViennaRNA 2004@ energy parameter file.

module Biobase.Turner.Import.Vienna where

import           Control.Applicative ( (<|>) )
import           Control.Monad (guard)
import           Data.List (foldl')
import           Data.Map.Strict (Map)
import           Data.PrimitiveArray as PA hiding (map,fromList)
import           Data.Semigroup ((<>))
import           Data.Text (Text,pack)
import           Data.Vector.Unboxed (Vector,fromList)
import qualified Data.Map.Strict as M
import           Text.Printf
import           Text.Trifecta as TT

import           Biobase.Primary (Primary,RNA,primary)
import           Biobase.Primary.Nuc.RNA
import           Biobase.Secondary.Vienna
import           Biobase.Turner.Types
import           Biobase.Types.Energy



-- |
--
-- TODO Use @parseFromFileEx@ and provide a function to display errors on
-- the shell.

fromFile :: FilePath -> IO (Maybe (Vienna2004,Vienna2004))
fromFile = TT.parseFromFile pVienna

pVienna :: Parser (Vienna2004,Vienna2004)
pVienna = v2Header *> spaces *> (blocksToTurner <$> some (block <* whiteSpace)) <* v2End <* eof
  where v2Header = text "## RNAfold parameter file v2.0" <?> ".par Header"
        v2End = text "#END" <* newline

-- | Run through the parsed blocks and insert into @Vienna 2004@ model.

blocksToTurner :: [Block] -> (Vienna2004,Vienna2004)
blocksToTurner = foldl' go (emptyViennaModel,emptyViennaModel)
  where go (s,h) = \case
          BlockPP "stack" Entropy  arr -> ( s {_stack = arr} , h )
          BlockPP "stack" Enthalpy arr -> ( s , h {_stack = arr} )
          BlockPBB t Entropy arr
            | t == "mismatch_hairpin"     -> (s {_hairpinMM = arr}  , h)
            | t == "mismatch_interior"    -> (s {_iloopMM = arr}    , h)
            | t == "mismatch_interior_1n" -> (s {_iloop1xnMM = arr} , h)
            | t == "mismatch_interior_23" -> (s {_iloop2x3MM = arr} , h)
            | t == "mismatch_multi"       -> (s {_multiMM = arr}    , h)
            | t == "mismatch_exterior"    -> (s {_exteriorMM = arr} , h)
          BlockPBB t Enthalpy arr
            | t == "mismatch_hairpin"     -> (s , h {_hairpinMM = arr} )
            | t == "mismatch_interior"    -> (s , h {_iloopMM = arr}   )
            | t == "mismatch_interior_1n" -> (s , h {_iloop1xnMM = arr})
            | t == "mismatch_interior_23" -> (s , h {_iloop2x3MM = arr})
            | t == "mismatch_multi"       -> (s , h {_multiMM = arr}   )
            | t == "mismatch_exterior"    -> (s , h {_exteriorMM = arr})
          BlockPB t Entropy arr
            | t == "dangle5" -> (s {_dangle5 = arr}, h )
            | t == "dangle3" -> (s {_dangle3 = arr}, h )
          BlockPB t Enthalpy arr
            | t == "dangle5" -> (s , h {_dangle5 = arr})
            | t == "dangle3" -> (s , h {_dangle3 = arr})
          BlockPPBB t Entropy arr
            | t == "int11" -> (s {_iloop1x1 = arr}, h )
          BlockPPBB t Enthalpy arr
            | t == "int11" -> (s , h {_iloop1x1 = arr})
          BlockPPBBB t Entropy arr
            | t == "int21" -> (s {_iloop2x1 = arr}, h )
          BlockPPBBB t Enthalpy arr
            | t == "int21" -> (s , h {_iloop2x1 = arr})
          BlockPPBBBB t Entropy arr
            | t == "int22" -> (s {_iloop2x2 = arr}, h )
          BlockPPBBBB t Enthalpy arr
            | t == "int22" -> (s , h {_iloop2x2 = arr})
          BlockLinear t Entropy arr
            | t == "hairpin"  -> (s {_hairpinL = arr} , h)
            | t == "bulge"    -> (s {_bulgeL = arr} , h)
            | t == "interior" -> (s {_bulgeL = arr} , h)
          BlockLinear t Enthalpy arr
            | t == "hairpin"  -> (s , h {_hairpinL = arr})
            | t == "bulge"    -> (s , h {_bulgeL = arr})
            | t == "interior" -> (s , h {_bulgeL = arr})
          -- | @u@ / @uh@ are unpaired nucleotides in the multibranched
          -- loop, as entropy and enthalpy. @b@ / @bl@ are the base cost
          -- for having a multibranched loop. @l@ / @lh@ are the gains for
          -- each stem in the multibranched loop.
          -- TODO write me
          BlockML t u uh b bh l lh ->
            (s , h)
          -- | @n@ / @nh@ are the ninio correction, @mx@ the maximal ninio
          -- correction
          BlockNinio t n nh mx ->
            (s {_ninio = n, _maxNinio = mx} , h {_ninio = nh})
          -- | @d@ / @dh@ is the duplex initiation energy, @au@ / @auh@
          -- is terminal-AU
          BlockMisc t d dh au auh ->
            (s {_termAU = au , _intermolecularInit = d }
            ,h {_termAU = auh, _intermolecularInit = dh})
          BlockLoops t rs ->
            (s {_hairpinLookup = _hairpinLookup s <> M.map fst rs}
            ,h {_hairpinLookup = _hairpinLookup h <> M.map snd rs})
          unknown       -> error $ "blocksToTurner: unknown " ++ show unknown

-- | Blocks have different arrays. The block parser will switch into the
-- corresponding parser where necessary.

data Block
  = BlockPP     Text EE (Arr VPP)
  | BlockPBB    Text EE (Arr VPNN)
  | BlockPB     Text EE (Arr VPN)
  | BlockPPBB   Text EE (Arr VPPNN)
  | BlockPPBBB  Text EE (Arr VPPNNN)
  | BlockPPBBBB Text EE (Arr VPPNNNN)
  | BlockLinear Text EE (Vector DeltaDekaGibbs)
  | BlockML     Text DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs
  | BlockNinio  Text DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs
  | BlockMisc   Text DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs DeltaDekaGibbs
  | BlockLoops  Text    (Map (Primary RNA) (DeltaDekaGibbs,DeltaDekaGibbs))
  deriving (Show)

data EE = Entropy | Enthalpy
  deriving (Show)

type Arr i = Unboxed i DeltaDekaGibbs

block :: Parser Block
block =   blockPP <|> blockPBB <|> blockPB <|> blockPPBB <|> blockPPBBB <|> blockPPBBBB
      <|> blockLinear <|> try blockML <|> try blockNinio <|> try blockMisc <|> blockLoops

-- | Parses a block of two stacked pairs.

blockPP :: Parser Block
blockPP = wrapBlock BlockPP blockHeader convertPP
  where convertPP vs = do
          let l = length vs
          guard (l == 7*7) <?> "expected 7*7=49 entries for blockPP but got " ++ show l
          let rs = fromAssocs minVPP maxVPP (DekaG 999999)
                 . guardLength l
                 $ zip [ (Z:.a:.b) | a <- cgnsP, b <- cgnsP ] vs
          return rs
        blockHeader = ["stack"]

guardLength k xs
  | k == length xs = xs
  | otherwise      = error $ printf "length xs /= %d in xs = %s" k $ show xs

acgun = [A,C,G,U,N]
nacgu = [N,A,C,G,U]

-- TODO we might have to introduce @E@ ...
eacgu = [N,A,C,G,U]

blockPBB ∷ Parser Block
blockPBB = wrapBlock BlockPBB blockHeader convertPBB
  where convertPBB vs = do
          let l = length vs
          guard (l == 7*5*5) <?> "expected 7*5*5=175 entries for blockPBB but got " ++ show l
          let rs = fromAssocs minVPBB maxVPBB (DekaG 999999)
                 . guardLength l
                 $ zip [ Z:.p:.a:.b | p ← cgnsP -- the closing pair
                       , a ← eacgu, b ← eacgu ] vs
          return rs
        blockHeader = ["mismatch_hairpin", "mismatch_interior", "mismatch_interior_1n", "mismatch_interior_23", "mismatch_multi", "mismatch_exterior"]

blockPB :: Parser Block
blockPB = wrapBlock BlockPB blockHeader convertPB
  where convertPB vs = do
          let l = length vs
          guard (l == 7*5) <?> "7*5=35 entries for blockPB but got " ++ show l
          let rs = fromAssocs minVPB maxVPB (DekaG 999999)
                 . guardLength l
                 $ zip [ Z:.p:.a | p ← cgnsP, a ← nacgu ] vs
          return rs
        blockHeader = ["dangle5", "dangle3"]

blockPPBB :: Parser Block
blockPPBB = wrapBlock BlockPPBB blockHeader convertPPBB
  where convertPPBB vs = do
          let l = length vs
          guard (l == 7*7*5*5) <?> "7*7*5*5=1225 entries for blockPPBB but got " ++ show l
          let rs = fromAssocs minVPPBB maxVPPBB (DekaG 999999)
                 . guardLength l
                 $ zip [ Z:.p:.q:.a:.b | p ← cgnsP, q ← cgnsP
                       , a ← nacgu, b ← nacgu ] vs
          return rs
        blockHeader = ["int11"]

blockPPBBB :: Parser Block
blockPPBBB = wrapBlock BlockPPBBB blockHeader convertPPBBB
  where convertPPBBB vs = do
          let l = length vs
          guard (l == 7*7*5*5*5) <?> "7*7*5*5*5=6125 entries for blockPPBBB but got " ++ show l
          let rs = fromAssocs minVPPBBB maxVPPBBB (DekaG 999999)
                 . guardLength l
                 $ zip [ Z:.p:.q:.a:.b:.c | p ← cgnsP, q ← cgnsP
                       , a ← nacgu, b ← nacgu, c ← nacgu ] vs
          return rs
        blockHeader = ["int21"]

-- | Parses the big @int22@ lists. Excludes special characters and only
-- includes standard nucleotides / pairs.

blockPPBBBB :: Parser Block
blockPPBBBB = wrapBlock BlockPPBBBB blockHeader convertPPBBBB
  where convertPPBBBB vs = do
          let l = length vs
          guard (l == 6*6*4*4*4*4) <?> "6*6*4*4*4*4=30625 entries for blockPPBBBB but got " ++ show l
          let rs = fromAssocs minVPPBBBB maxVPPBBBB (DekaG 999999)
                 . guardLength l
                 $ zip [ Z:.p:.q:.a:.b:.c:.d | p ← cguaP, q ← cguaP
                       , a ← acgu, b ← acgu, c ← acgu, d ← acgu ] vs
          return rs
        blockHeader = ["int22"]

blockLinear :: Parser Block
blockLinear = wrapBlock BlockLinear blockHeader convertLinear
  where convertLinear vs = do
          let l = length vs
          guard (l == 31) <?> "31 entries for blockLinear but got " ++ show l
          return $ fromList vs
        blockHeader = ["hairpin", "bulge", "interior"]

blockML :: Parser Block
blockML = BlockML <$  text "# " <*> text "ML_params" <* newline
                  <*> ddg <*> ddg <*> ddg <*> ddg <*> ddg <*> ddg

blockNinio :: Parser Block
blockNinio = BlockNinio <$  text "# " <*> text "NINIO" <* newline
                        <*> ddg <*> ddg <*> ddg

blockMisc :: Parser Block
blockMisc = BlockMisc <$  text "# " <*> text "Misc" <* newline
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
  where pEE = (Entropy <$ newline) <|> (Enthalpy <$ text "_enthalpies" <* newline)

arrayLines :: Parser [DeltaDekaGibbs]
arrayLines = do
  headerComment <- commentLine
  valueLines <- valueLine `sepEndBy` newline
  return . map DekaG $ concat valueLines

commentLine = optional $ between (text "/*") (text "*/") (spaces *> some (letter <|> oneOf ",") `sepEndBy` spaces) <?> "commentLine"

valueLine :: Parser [Int]
valueLine = spaces *> (fromInteger <$> integer <|> 999999 <$ text "INF") `sepEndBy` spaces <* commentLine -- <* newline


test = fromFile "data/rna_turner2004.par"

