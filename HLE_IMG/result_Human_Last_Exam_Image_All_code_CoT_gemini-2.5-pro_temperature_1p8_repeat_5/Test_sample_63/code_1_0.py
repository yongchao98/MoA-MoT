def solve_nmr():
    """
    Analyzes the provided 1H NMR spectrum and determines the corresponding molecular structure
    from the given candidates.
    """
    print("--- 1H NMR Spectrum Analysis ---")
    print("\nStep 1: Signal-by-signal breakdown of the spectrum.")

    # Signal 1: ~1.1 ppm
    print("\nSignal 1 (at ~1.1 ppm):")
    print("  - Chemical Shift: ~1.1 ppm (upfield, alkane region)")
    print("  - Multiplicity: Triplet (split by 2 neighboring protons, -CH2-)")
    print("  - Integration: Corresponds to 6 protons (relative ratio).")
    print("  - Deduced Fragment: Two equivalent methyl groups adjacent to a CH2 group, i.e., -N(CH2-CH3)2.")

    # Signal 2: ~2.2 ppm
    print("\nSignal 2 (at ~2.2 ppm):")
    print("  - Chemical Shift: ~2.2 ppm (benzylic region)")
    print("  - Multiplicity: Singlet (no neighboring protons)")
    print("  - Integration: Corresponds to 3 protons.")
    print("  - Deduced Fragment: An isolated methyl group attached to the aromatic ring (Ar-CH3).")

    # Signal 3: ~2.7 ppm
    print("\nSignal 3 (at ~2.7 ppm):")
    print("  - Chemical Shift: ~2.7 ppm")
    print("  - Multiplicity: Quartet (split by 3 neighboring protons, -CH3)")
    print("  - Integration: Corresponds to 4 protons.")
    print("  - Deduced Fragment: Two equivalent methylene groups adjacent to a CH3 group. This confirms the diethylamino fragment, -N(CH2-CH3)2.")

    # Signal 4: ~3.3 ppm
    print("\nSignal 4 (at ~3.3 ppm):")
    print("  - Chemical Shift: ~3.3 ppm (downfield, next to electronegative atoms)")
    print("  - Multiplicity: Singlet (no neighboring protons)")
    print("  - Integration: Corresponds to 2 protons.")
    print("  - Deduced Fragment: An isolated methylene group between a carbonyl and a nitrogen, -C(=O)-CH2-N-.")

    # Signal 5: ~7.0-7.5 ppm
    print("\nSignal 5 (at ~7.0-7.5 ppm):")
    print("  - Chemical Shift: ~7.2 ppm (aromatic region)")
    print("  - Multiplicity: Multiplet (complex splitting from other ring protons)")
    print("  - Integration: Corresponds to 3 protons.")
    print("  - Deduced Fragment: Three protons on a substituted benzene ring.")

    # Signal 6: ~8.9 ppm
    print("\nSignal 6 (at ~8.9 ppm):")
    print("  - Chemical Shift: ~8.9 ppm (very downfield)")
    print("  - Multiplicity: Broad Singlet")
    print("  - Integration: Corresponds to 1 proton.")
    print("  - Deduced Fragment: An amide proton, -NH-C(=O)-.")

    print("\n--- Step 2: Evaluating Candidate Structures ---")

    print("\nCandidate A-G and B-G:")
    print("  - These have a dimethylamino group, -N(CH3)2, which would give a 6H singlet.")
    print("  - The spectrum shows a 6H triplet and 4H quartet, not a 6H singlet.")
    print("  - Conclusion: A-G and B-G are incorrect.")

    print("\nCandidate D-L:")
    print("  - This structure has two methyl groups on the aromatic ring, which would give signals for 6 protons (e.g., a 6H singlet or two 3H singlets) around 2.2 ppm.")
    print("  - The spectrum shows only one 3H singlet in this region.")
    print("  - Conclusion: D-L is incorrect.")

    print("\nCandidate C-L:")
    print("  - This structure's features match the spectrum perfectly:")
    print(f"  - -N(C2H5)2 group: Explains the 6H triplet (~1.1 ppm) and 4H quartet (~2.7 ppm).")
    print(f"  - Ar-CH3 group: Explains the 3H singlet at ~2.2 ppm.")
    print(f"  - -C(=O)-CH2-N- group: Explains the 2H singlet at ~3.3 ppm.")
    print(f"  - Aromatic ring protons: A trisubstituted ring has 3H, matching the multiplet at ~7.2 ppm.")
    print(f"  - Amide -NH proton: Explains the 1H broad singlet at ~8.9 ppm.")
    print("  - Conclusion: C-L is the correct structure.")

solve_nmr()