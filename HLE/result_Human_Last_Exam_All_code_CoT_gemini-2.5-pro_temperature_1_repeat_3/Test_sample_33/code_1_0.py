def analyze_protein_folding():
    """
    Analyzes FTIR data to explain the structural changes in tardigrade proteins upon gelation.
    """

    # Step 1: Assign FTIR peaks to their corresponding secondary structures.
    print("Step 1: Assigning FTIR peaks to protein secondary structures.")
    peak_assignments = {
        1645: "Disordered / Random Coil (typically broad)",
        1652: "Alpha-Helix",
        1618: "Beta-Sheet",
        1680: "Beta-Sheet (often indicating anti-parallel beta-sheets)"
    }
    for peak, structure in peak_assignments.items():
        print(f" - A peak at {peak} cm^-1 corresponds to: {structure}")
    print("-" * 30)

    # Step 2: Analyze the concentration titration experiment, which shows the gelation process.
    print("Step 2: Analyzing the concentration titration (gelation).")
    print("The initial state of the protein is described as 'disordered', consistent with the presence of a peak around 1645 cm^-1.")
    print("As concentration increases, gelation occurs.")
    print(f" - The observation of an increase in the {1652} cm^-1 peak indicates the FORMATION of alpha-helices.")
    print(f" - The simultaneous increase in the {1618} cm^-1 peak indicates the FORMATION of beta-sheets.")
    print("Conclusion from this experiment: Gelation involves disordered structures folding into both alpha-helices and beta-sheets.")
    print("-" * 30)

    # Step 3: Use the heating experiment to confirm the findings. Heating reverses gelation.
    print("Step 3: Analyzing the heating experiment (de-gelation).")
    print("Upon heating, the gel is expected to denature or unfold back towards a disordered state.")
    print(f" - The disappearance of the {1618} cm^-1 and {1680} cm^-1 peaks confirms that beta-sheets are part of the ordered gel structure that is lost upon heating.")
    print(f" - The strengthening of the {1645} cm^-1 peak confirms the transition to a more disordered state.")
    print("This confirms that the gel state is an ordered structure containing beta-sheets, which unfolds into a disordered structure when heated.")
    print("-" * 30)

    # Step 4: Synthesize and select the final answer.
    print("Step 4: Final Conclusion.")
    print("Both experiments support the same conclusion: The initially disordered proteins fold to form a gel structure that contains both alpha-helical and beta-sheet content.")
    print("This matches the answer choice describing the folding of disordered structures into both beta sheets and alpha helices.")

analyze_protein_folding()