def find_molecule_from_nmr():
    """
    Analyzes a 1H NMR spectrum to identify the correct molecular structure
    from a list of candidates by systematically matching spectral features
    to structural fragments.
    """

    print("Step 1: Analysis of the 1H NMR Spectrum")
    print("The spectrum shows the following signals:")
    print(" - A triplet at ~1.1 ppm (Integration: 6H)")
    print(" - A singlet at ~2.3 ppm (Integration: 3H)")
    print(" - A quartet at ~2.8 ppm (Integration: 4H)")
    print(" - A singlet at ~3.5 ppm (Integration: 2H)")
    print(" - A multiplet in the aromatic region at ~7.1 ppm (Integration: ~4H)")
    print(" - A broad singlet at ~8.9 ppm (Integration: 1H)")
    print("\n")

    print("Step 2: Deducing Structural Fragments from the Spectrum")
    print("The combination of the triplet at 1.1 ppm (6H) and the quartet at 2.8 ppm (4H) is a classic pattern for two ethyl groups attached to a heteroatom, such as a diethylamino group [-N(CH2CH3)2].")
    print("This observation immediately eliminates structures A-G and B-G, as they both contain a dimethylamino group [-N(CH3)2], which would show a single 6H singlet, not a triplet and a quartet.")
    print("\n")

    print("Step 3: Comparing the remaining candidates (C-L and D-L)")
    print("Structure D-L has two different methyl groups on the benzene ring. This would result in two separate 3H singlets in the NMR spectrum. However, we only observe one 3H singlet at 2.3 ppm. Therefore, structure D-L is inconsistent with the data.")
    print("\n")

    print("Step 4: Confirming the match with structure C-L")
    print("Structure C-L is fully consistent with all observed signals. Here is the final assignment:")
    print("===============================================================")
    print("Signal (ppm) | Protons | Multiplicity | Assignment in C-L")
    print("---------------------------------------------------------------")
    print(f"      8.9      |   1H    |    Singlet   | Amide N-H")
    print(f"      7.1      |   4H    |   Multiplet  | Aromatic H's on the tolyl ring")
    print(f"      3.5      |   2H    |    Singlet   | -C(=O)-CH2-N group")
    print(f"      2.8      |   4H    |    Quartet   | -N-(CH2-CH3)2 group")
    print(f"      2.3      |   3H    |    Singlet   | -CH3 group on the aromatic ring")
    print(f"      1.1      |   6H    |    Triplet   | -N-(CH2-CH3)2 group")
    print("===============================================================")
    print("\nConclusion: The NMR spectrum perfectly matches the structure of C-L.")


find_molecule_from_nmr()