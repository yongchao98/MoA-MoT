def analyze_ftir_data():
    """
    Analyzes FTIR data of tardigrade hydrogel proteins and explains the observations.
    """
    
    # Step 1: Define the FTIR peaks and their structural assignments
    peak_assignments = {
        1645: "Random Coil / Disordered Structure (Broad Peak)",
        1652: "Alpha Helix",
        1618: "Intermolecular Beta Sheet",
        1680: "Antiparallel Beta Sheet (often paired with a ~1620-1640 cm^-1 peak)"
    }
    
    print("--- FTIR Peak Analysis for Tardigrade Protein Gelation ---")
    print("\nStep 1: Assigning spectral peaks to protein secondary structures.")
    for peak, structure in peak_assignments.items():
        print(f"  - The peak at {peak} cm^-1 corresponds to: {structure}")

    print("\nStep 2: Interpreting the experimental observations.")
    
    print("\n  Observation A: Heating Experiment")
    print(f"  - The peak at 1645 cm^-1 (Disordered) grows stronger.")
    print(f"  - The peaks at 1618 cm^-1 (Beta Sheet) and 1680 cm^-1 (Beta Sheet) disappear.")
    print("  - Interpretation: This indicates thermal denaturation. The ordered beta-sheet structures are disrupted by heat and unfold into disordered random coils.")

    print("\n  Observation B: Concentration Titration Experiment")
    print(f"  - A dual increase is seen in the peaks at 1652 cm^-1 (Alpha Helix) and 1618 cm^-1 (Beta Sheet) as concentration increases.")
    print("  - Interpretation: This experiment mimics the gelation process. As proteins become more concentrated, they self-assemble.")
    
    print("\nStep 3: Synthesizing the information to explain gelation.")
    print("  - The initial state of the protein is described as disordered, which aligns with the strong 1645 cm^-1 signal.")
    print("  - The concentration experiment shows that upon gelation, the proteins don't just form one type of structure.")
    print(f"  - Instead, there is a clear and simultaneous increase in signals for BOTH Alpha Helices ({1652} cm^-1) and Beta Sheets ({1618} cm^-1).")

    print("\n--- Conclusion ---")
    print("The most accurate explanation is that the initially disordered proteins fold into a hydrogel network composed of a mix of both alpha-helical and beta-sheet structures.")
    print("This matches Answer Choice I.")

# Execute the analysis
analyze_ftir_data()