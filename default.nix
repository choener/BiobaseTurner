{ mkDerivation, aeson, base, binary, BiobaseTypes, BiobaseXNA
, bytestring, bytestring-lexing, cereal, containers, data-default
, file-embed, filepath, fused-effects, lens, lib, parsers
, primitive, PrimitiveArray, SciBaseTypes, split, text, trifecta
, vector, vector-th-unbox
}:
mkDerivation {
  pname = "BiobaseTurner";
  version = "0.3.2.1";
  src = ./.;
  libraryHaskellDepends = [
    aeson base binary BiobaseTypes BiobaseXNA bytestring
    bytestring-lexing cereal containers data-default file-embed
    filepath fused-effects lens parsers primitive PrimitiveArray
    SciBaseTypes split text trifecta vector vector-th-unbox
  ];
  homepage = "https://github.com/choener/BiobaseTurner";
  description = "The Turner RNA energy model";
  license = lib.licenses.gpl3Plus;
}
