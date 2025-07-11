def solve_nmr_puzzle():
    """
    This function provides a step-by-step logical deduction to identify the correct
    molecular structure based on the provided 1H NMR spectrum.
    """

    print("Step 1: Analyzing the 1H NMR Spectrum Features")
    print("-------------------------------------------------")
    print("Based on the spectrum, we can identify the following signals:")
    print("- A broad singlet at ~9.0 ppm (likely an N-H proton).")
    print("- A multiplet in the aromatic region at ~7.2 ppm.")
    print("- A sharp singlet at ~3.4 ppm (a CH2 group with no H neighbors).")
    print("- A singlet at ~2.4 ppm (a CH3 group with no H neighbors, like Ar-CH3).")
    print("- A quartet at ~2.8 ppm and a triplet at ~1.2 ppm.")
    print("\nThe quartet and triplet signals are a clear indicator of an ethyl group (-CH2CH3).\n")

    print("Step 2: Evaluating the Candidate Structures")
    print("------------------------------------------")
    print("A) Eliminating Candidates A-G and B-G:")
    print("   - Both A-G and B-G have a dimethylamino [-N(CH3)2] group, not a diethylamino group.")
    print("   - They would show a large 6H singlet for the two methyls, not the quartet and triplet pattern seen in the spectrum.")
    print("   - Therefore, candidates A-G and B-G are incorrect.\n")
    
    print("B) Differentiating Between C-L and D-L:")
    print("   - Both C-L and D-L have the diethylamino [-N(CH2CH3)2] group, matching the quartet at 2.8 ppm and triplet at 1.2 ppm.")
    print("   - The key difference is the aromatic ring.")
    print("   - Structure C-L: Has ONE aromatic methyl group (3 protons) and FOUR aromatic protons.")
    print("     Expected proton ratio for Ar-H (~7.2 ppm) vs. Ar-CH3 (~2.4 ppm) is 4 : 3.")
    print("   - Structure D-L: Has TWO aromatic methyl groups (6 protons) and THREE aromatic protons.")
    print("     Expected proton ratio for Ar-H (~7.2 ppm) vs. Ar-CH3 (~2.4 ppm) is 3 : 6 (or 1 : 2).\n")

    print("Step 3: Final Conclusion")
    print("------------------------")
    print("By visually inspecting the spectrum:")
    print("- The signal for the aromatic protons at 7.2 ppm is visibly larger than the signal for the aromatic methyl group at 2.4 ppm.")
    print("- This observation matches the predicted 4:3 proton ratio of structure C-L.")
    print("- It directly contradicts the 1:2 ratio expected for structure D-L.")
    print("\nTherefore, the spectrum corresponds to structure C-L.")

solve_nmr_puzzle()