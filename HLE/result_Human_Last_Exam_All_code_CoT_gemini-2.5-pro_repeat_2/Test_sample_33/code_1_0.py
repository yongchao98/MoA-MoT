def analyze_protein_folding():
    """
    Analyzes FTIR data to determine protein structural changes during gelation.
    """
    # Step 1: Define FTIR peak assignments based on known literature values.
    peak_assignments = {
        1652: "Alpha Helix",
        1618: "Beta Sheet",
        1645: "Disordered / Random Coil",
        1680: "Beta Sheet (anti-parallel component)"
    }

    print("--- Analysis of Tardigrade Protein Gelation ---")
    print("\nInitial State: Proteins are described as 'initially disordered'.")
    print("This corresponds to the structure typically seen around 1645 cm^-1.")

    print("\nObservation during Concentration Titration (Gelation):")
    print("A dual increase is observed for two specific peaks.")

    # Step 2: Analyze the observed changes during gelation.
    gelation_peak_1 = 1652
    gelation_peak_2 = 1618

    structure_formed_1 = peak_assignments[gelation_peak_1]
    structure_formed_2 = peak_assignments[gelation_peak_2]

    print(f"\n1. The peak at {gelation_peak_1} cm^-1 increases.")
    print(f"   - This signal corresponds to the formation of '{structure_formed_1}' structures.")

    print(f"\n2. The peak at {gelation_peak_2} cm^-1 also increases.")
    print(f"   - This signal corresponds to the formation of '{structure_formed_2}' structures.")

    # Step 3: Synthesize the information to form a conclusion.
    print("\n--- Conclusion ---")
    print("The experiment shows that upon increasing concentration, the initially disordered proteins")
    print(f"fold into a mixture of secondary structures. The evidence is:")
    print(f"  - Increase at {gelation_peak_1} cm^-1 => Formation of {structure_formed_1}")
    print(f"  - Increase at {gelation_peak_2} cm^-1 => Formation of {structure_formed_2}")
    print("\nTherefore, the most accurate description is that disordered structures fold into both beta sheets and alpha helices upon gelation.")

if __name__ == '__main__':
    analyze_protein_folding()
<<<I>>>