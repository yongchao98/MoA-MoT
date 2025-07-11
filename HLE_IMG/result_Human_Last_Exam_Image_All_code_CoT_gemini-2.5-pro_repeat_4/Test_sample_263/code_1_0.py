def analyze_lec_stability():
    """
    Analyzes the expected stability of LECs based on three Ir(III) complexes
    by assigning a score based on known molecular design principles.
    """

    # Define stability features for each complex:
    # [has_fluorination, has_steric_bulk, has_reactive_ligand]
    complex_features = {
        "Complex 1": [False, False, False],
        "Complex 2": [False, False, True],
        "Complex 3": [True, True, False]
    }

    stability_scores = {}
    print("Analyzing the expected device stability for each complex...")
    print("="*60)
    print("Scoring Rules:")
    print(" - Base Score: 1")
    print(" - Fluorination on C^N ligand: +1 (enhances bond strength/oxidative stability)")
    print(" - Steric Bulk on N^N ligand: +1 (protects from intermolecular degradation)")
    print(" - Large, electroactive N^N ligand: -1 (potential degradation pathway)")
    print("="*60)

    for name, features in complex_features.items():
        # Each complex starts with a base stability score of 1
        score = 1
        has_fluorination, has_steric_bulk, has_reactive_ligand = features

        if has_fluorination:
            score += 1
        if has_steric_bulk:
            score += 1
        if has_reactive_ligand:
            score -= 1
        
        stability_scores[name] = score

    # Print the analysis and scores
    print("\nStability Analysis and Scores:")
    print(f"Complex 1 (Baseline): Score = 1")
    print(f"Complex 2 (Reactive Ligand): Score = 1 - 1 = 0")
    print(f"Complex 3 (Fluorination + Steric Bulk): Score = 1 + 1 + 1 = 3")
    print("-" * 35)

    # Determine the most stable complex
    most_stable_complex = max(stability_scores, key=stability_scores.get)

    print("\nConclusion:")
    print(f"Based on the analysis, {most_stable_complex} has the highest stability score.")
    print("The combination of fluorinated cyclometalating ligands and a sterically hindered ancillary ligand in Complex 3")
    print("provides superior protection against common degradation mechanisms.")
    print("\nTherefore, LECs based on Complex 3 are expected to result in more stable devices.")

if __name__ == '__main__':
    analyze_lec_stability()