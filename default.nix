{ mkDerivation, aeson, base, binary, BiobaseTypes, BiobaseXNA
, bytestring, bytestring-lexing, cereal, cmdargs, containers
, data-default, file-embed, filepath, fused-effects, HUnit, lens
, primitive, PrimitiveArray, QuickCheck, SciBaseTypes, split
, stdenv, tasty, tasty-hunit, tasty-quickcheck, tasty-th, text
, trifecta, vector, vector-th-unbox
}:
mkDerivation {
  pname = "BiobaseTurner";
  version = "0.3.2.1";
  src = ./.;
  isLibrary = true;
  isExecutable = true;
  libraryHaskellDepends = [
    aeson base binary BiobaseTypes BiobaseXNA bytestring
    bytestring-lexing cereal containers data-default file-embed
    filepath fused-effects lens primitive PrimitiveArray SciBaseTypes
    split text trifecta vector vector-th-unbox
  ];
  executableHaskellDepends = [
    aeson base binary bytestring cereal cmdargs containers vector
  ];
  testHaskellDepends = [
    base BiobaseXNA data-default HUnit lens PrimitiveArray QuickCheck
    tasty tasty-hunit tasty-quickcheck tasty-th
  ];
  homepage = "https://github.com/choener/BiobaseTurner";
  description = "Import Turner RNA parameters";
  license = stdenv.lib.licenses.gpl3;
}
