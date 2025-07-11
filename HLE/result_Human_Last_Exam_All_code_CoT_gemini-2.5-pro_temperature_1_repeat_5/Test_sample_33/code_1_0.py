def analyze_protein_folding():
    """
    Analyzes FTIR data to explain the structural changes of a tardigrade protein upon hydrogel formation.
    """

    # Step 1: Define the standard FTIR peak assignments for protein secondary structures.
    print("Step 1: Assigning FTIR peaks to protein secondary structures based on their wavenumbers.")
    print("--------------------------------------------------------------------------------------")
    print(f"Peak at 1652 cm^-1: Corresponds to Alpha-Helices.")
    print(f"Peak at 1618 cm^-1: Corresponds to Beta-Sheets.")
    print(f"Peak at 1680 cm^-1: Corresponds to anti-parallel Beta-Sheets (a specific type).")
    print(f"Peak at 1645 cm^-1 (broad): Corresponds to Disordered/Random Coil structures.")
    print("\n")

    # Step 2: Analyze the concentration titration experiment, which shows how the gel forms.
    print("Step 2: Analyzing the concentration titration experiment.")
    print("-------------------------------------------------------")
    print("Observation: As protein concentration increases, there is a dual increase in the peaks at 1652 cm^-1 and 1618 cm^-1.")
    print("Interpretation: The protein starts as disordered. The increase in these two peaks means that upon gelation, the disordered protein chains are folding into both Alpha-Helices (from the 1652 cm^-1 peak) and Beta-Sheets (from the 1618 cm^-1 peak).")
    print("\n")

    # Step 3: Use the heating experiment to confirm the nature of the structures.
    print("Step 3: Analyzing the heating experiment for confirmation.")
    print("---------------------------------------------------------")
    print("Observation: Upon heating, the Beta-Sheet peaks (1618 cm^-1 and 1680 cm^-1) disappear, while the Disordered peak (1645 cm^-1) grows stronger.")
    print("Interpretation: This confirms that the Beta-Sheets are part of the ordered, heat-sensitive gel structure which unfolds back into a disordered state with heat.")
    print("\n")

    # Step 4: Synthesize the findings to reach a final conclusion.
    print("Step 4: Final Conclusion.")
    print("-------------------------")
    print("The concentration experiment is key: it shows that the gel is FORMED by creating both alpha-helices and beta-sheets from an initially disordered state.")
    print("Therefore, the most accurate explanation for the behavior is that the disordered protein structure folds into both beta sheets and alpha helices upon gelation.")

# Execute the analysis function
analyze_protein_folding()
