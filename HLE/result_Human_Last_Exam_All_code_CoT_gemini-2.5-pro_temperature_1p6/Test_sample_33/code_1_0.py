import textwrap

def analyze_protein_folding():
    """
    Analyzes FTIR data to determine the protein folding behavior upon gelation.
    """

    # Step 1: Define FTIR peak assignments
    peak_assignments = {
        1645: "Disordered / Random Coil (Broad peak)",
        1652: "Alpha Helix",
        1618: "Intermolecular Beta Sheet (Low frequency)",
        1680: "Antiparallel Beta Sheet (High frequency component)"
    }

    # Helper function for pretty printing
    def print_wrapped(text):
        print(textwrap.fill(text, width=80))

    print("### Analysis of Tardigrade Protein Gelation ###\n")

    # Step 2: Print peak assignments for clarity
    print_wrapped("1. Assigning Protein Structures to FTIR Peaks:")
    for peak, structure in peak_assignments.items():
        print(f"   - {peak} cm⁻¹: Corresponds to {structure}")
    print("-" * 20)

    # Step 3: Analyze the experimental observations
    print_wrapped("2. Analyzing the Experimental Data:")

    print("\n--- Initial State ---")
    print_wrapped(f"The protein is initially a mix of structures. The strong, broad peak at {1645} cm⁻¹ indicates a large portion of the protein is in a Disordered state. The shoulder peaks at {1652}, {1618}, and {1680} cm⁻¹ indicate the presence of some pre-existing Alpha Helices and Beta Sheets.")

    print("\n--- Heating Experiment ---")
    print_wrapped(f"Observation: Upon heating, the {1645} cm⁻¹ (Disordered) peak grows, while the {1618} and {1680} cm⁻¹ (Beta Sheet) peaks disappear.")
    print_wrapped(f"Interpretation: Heat causes denaturation. The ordered structures (specifically Beta Sheets) are unfolding into a Disordered state. This is a typical behavior for proteins.")

    print("\n--- Concentration Titration (Gelation) ---")
    print_wrapped(f"Observation: As concentration increases and the hydrogel forms, there is a dual increase in the {1652} cm⁻¹ (Alpha Helix) peak and the {1618} cm⁻¹ (Beta Sheet) peak.")
    print_wrapped("Interpretation: Gelation is a process of ordering. Since the protein is initially disordered, this observation shows that the disordered protein chains are folding into BOTH Alpha Helices AND Beta Sheets to form the gel network.")
    print("-" * 20)

    # Step 4: Evaluate the options
    print("3. Evaluating the Answer Choices:\n")
    print_wrapped("A, D, E, H are incorrect because they describe unfolding or transformation from one ordered state to another, but the data clearly shows a net INCREASE in ordered structures (both alpha-helix and beta-sheet) from a disordered state upon gelation.")
    print_wrapped("B and F are incomplete. They correctly identify folding from a disordered state but only mention one of the two structures observed to be forming.")
    print_wrapped("G is incorrect as it describes a rearrangement of existing beta-sheets, not the formation of both helices and sheets from a disordered state.")
    print_wrapped("I is the most accurate description. It accounts for the initial disordered state and the observed, simultaneous increase of both alpha-helical and beta-sheet structures during concentration-induced gelation.")

    # Step 5: Conclude
    print("\n### Conclusion ###")
    print_wrapped("The evidence strongly supports that the gelation process involves the protein folding from an initially disordered state into a final structure containing both alpha-helices and beta-sheets.")

analyze_protein_folding()
