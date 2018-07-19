with (import <nixpkgs> {});
with haskell.lib;

rec {
  hsSrcSet = (lib.foldl' (s: p: s // (import p).hsSrcSet) {} [
    ../Lib-BiobaseTypes
    ../Lib-BiobaseXNA
    ../Lib-PrimitiveArray
    ../Lib-SciBaseTypes
  ]) // {BiobaseTurner = ./.;};
  hsPkgs = haskellPackages.extend (packageSourceOverrides hsSrcSet);
  hsShell = with hsPkgs; shellFor {
    packages = p: [ p.BiobaseTurner ];
    withHoogle = true;
    buildInputs = [
      cabal-install ghc
      BiobaseTypes
      BiobaseXNA
      PrimitiveArray
      SciBaseTypes
    ];
  };
}
