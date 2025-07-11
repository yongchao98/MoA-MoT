def solve_nmr_puzzle():
    """
    Analyzes the provided 13C NMR data for a C7H14 hydrocarbon to determine its IUPAC name.
    """
    # --- Step 1: Define the given data ---
    molecular_formula = "C7H14"
    nmr_signals = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }
    num_carbons_formula = 7
    num_signals = len(nmr_signals)

    # --- Step 2: Print the analysis ---
    print("--- Analysis of the Spectroscopic Data ---")
    print(f"Molecular Formula: {molecular_formula}")

    # Calculate Degree of Unsaturation (DoU)
    h_saturated = 2 * num_carbons_formula + 2
    h_actual = 14
    dou = (h_saturated - h_actual) // 2
    print(f"Degree of Unsaturation: {dou}. This indicates one double bond or one ring.")

    # Analyze signal count vs carbon count
    print(f"Number of 13C NMR signals: {num_signals}")
    print(f"Number of Carbon atoms: {num_carbons_formula}")
    print("Since there are 6 signals for 7 carbons, the molecule has symmetry, and two carbons are chemically equivalent.")
    print("\n--- Interpretation of Individual Signals ---")

    # Analyze each signal
    print("Signal at 145 ppm (s): In alkene region (100-150 ppm), 's' indicates a quaternary carbon (0H). This is a >C= carbon.")
    print("Signal at 112 ppm (t): In alkene region, 't' indicates a CH2 group (2H). This is a terminal =CH2 group.")
    print("Conclusion from alkene signals: The molecule contains a >C=CH2 fragment, confirming a double bond.")
    print("\nSignal at 48 ppm (t): Aliphatic region, 't' indicates a CH2 group.")
    print("Signal at 27 ppm (d): Aliphatic region, 'd' indicates a CH group.")
    print("Signal at 22 ppm (q): Aliphatic region, 'q' indicates a CH3 group.")
    print("Signal at 21 ppm (q): Aliphatic region, 'q' indicates another, different CH3 group.")

    # --- Step 3: Propose a structure and give the final answer ---
    print("\n--- Structure Determination ---")
    print("The fragments (1x >C=, 1x =CH2, 1x -CH2-, 1x -CH-, and two different -CH3 environments) can be assembled into a single structure.")
    print("The structure consistent with all data is 2-ethyl-3-methyl-1-butene.")
    print("This structure has 7 carbons and 14 hydrogens, matching the formula.")
    print("It has 6 unique carbon environments, which matches the number of signals:")
    print(" - C1 (=CH2): Corresponds to 112(t)")
    print(" - C2 (>C=): Corresponds to 145(s)")
    print(" - Ethyl -CH2-: Corresponds to 48(t)")
    print(" - C3 (-CH-): Corresponds to 27(d)")
    print(" - Ethyl -CH3: Corresponds to one of the signals 22(q) or 21(q)")
    print(" - The two equivalent methyls on C3: Corresponds to the other 'q' signal, 22(q) or 21(q)")

    print("\n--- Final Answer ---")
    final_iupac_name = "2-ethyl-3-methyl-1-butene"
    print(f"The IUPAC name of the compound is: {final_iupac_name}")

if __name__ == '__main__':
    solve_nmr_puzzle()