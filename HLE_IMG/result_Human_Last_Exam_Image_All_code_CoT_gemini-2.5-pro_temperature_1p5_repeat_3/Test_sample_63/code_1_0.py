def solve_nmr_puzzle():
    """
    This program analyzes the provided 1H NMR data to determine the correct structure
    from the given candidates by assigning peaks to proton environments.
    """

    print("Step 1: Initial analysis of the NMR spectrum.")
    print("The spectrum shows a triplet at ~1.0 ppm and a quartet at ~2.5 ppm.")
    print("This pattern is a classic signature for ethyl groups (-CH2-CH3).")
    print("The relative integrations (quartet:triplet ratio of 4H:6H) confirm the presence of a diethylamino group: -N(CH2-CH3)2.")
    print("This key feature eliminates candidates A-G and B-G, which contain dimethylamino (-N(CH3)2) groups.")

    print("\nStep 2: Differentiating between the remaining candidates, C-L and D-L.")
    print("Both C-L and D-L contain the diethylamino group and a methyl-substituted benzene ring with an amide linker.")
    print("We must look for differences in their expected spectra.")

    print("\nPredicted signals for both C-L and D-L:")
    # Assignments are based on standard chemical shift values and splitting rules.
    peak_assignments = {
        "N(CH2-CH3)2 (methyl protons, 6H)": "A triplet at approx. 1.0 ppm. Matches spectrum.",
        "Ar-CH3 (aromatic methyl, 3H)": "A singlet at approx. 2.2 ppm. Matches spectrum.",
        "N(CH2-CH3)2 (methylene protons, 4H)": "A quartet at approx. 2.5 ppm. Matches spectrum.",
        "-CO-CH2-N- (methylene protons, 2H)": "A singlet at approx. 3.4 ppm. Matches spectrum.",
        "Aromatic Protons (4H for C-L, 4H for D-L)": "Signals around 7.1 ppm. Matches spectrum.",
        "Amide N-H (1H)": "A singlet at approx. 8.8 ppm. Matches spectrum."
    }
    for proton_group, description in peak_assignments.items():
      print(f"- {proton_group}: {description}")


    print("\nStep 3: The deciding factor - the aromatic region (~7.1 ppm).")
    print("Structure D-L is para-substituted (methyl at C4). Due to symmetry, it would show a simpler pattern for its 4 aromatic protons (an AA'BB' system, often appearing as two doublets).")
    print("Structure C-L is ortho-substituted (methyl at C2). It has no such symmetry, so its 4 aromatic protons are all different, giving a complex, overlapping multiplet (an ABCD system).")
    print("The observed spectrum shows a complex multiplet at ~7.1 ppm, which is consistent with the ortho-substituted structure C-L.")

    print("\n--- CONCLUSION ---")
    print("The most possible structure is C-L.")
    print("The full assignment for C-L matches the spectrum:")
    print("triplet at 1.0 ppm (6H) + quartet at 2.5 ppm (4H) = -N(C2H5)2 group")
    print("singlet at 2.2 ppm (3H) = Ar-CH3 group")
    print("singlet at 3.4 ppm (2H) = -CO-CH2-N- group")
    print("complex multiplet at 7.1 ppm (4H) = ortho-substituted aromatic ring")
    print("singlet at 8.8 ppm (1H) = -NH- amide proton")


solve_nmr_puzzle()