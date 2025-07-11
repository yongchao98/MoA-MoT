def analyze_emitter_stability():
    """
    Analyzes the expected stability of three Ir(III) complexes in LECs
    based on their structural features.
    """

    # 1. Represent complexes with their key structural features
    complexes = {
        'Complex 1': {
            'cyclometalating_ligand': 'ppy (unsubstituted)',
            'ancillary_ligand': 'bpy (unsubstituted)',
            'features': []
        },
        'Complex 2': {
            'cyclometalating_ligand': 'ppy (unsubstituted)',
            'ancillary_ligand': 'Bulky phenanthro-imidazole derivative',
            'features': ['bulky_groups']
        },
        'Complex 3': {
            'cyclometalating_ligand': 'dfppy (fluorinated)',
            'ancillary_ligand': 'dtbbpy (bulky)',
            'features': ['fluorination', 'bulky_groups']
        }
    }

    # 2. Define stability scores based on chemical principles
    stability_scores = {
        'fluorination': 2,
        'bulky_groups': 1
    }

    print("Evaluating LEC stability based on emitter molecular design:")
    print("-" * 60)
    print("Scoring Rules:")
    print(f"- Fluorination of cyclometalating ligand: +{stability_scores['fluorination']} points (enhances bond strength)")
    print(f"- Bulky groups on ligands: +{stability_scores['bulky_groups']} point (provides steric protection)\n")

    final_scores = {}

    # 3. Calculate and print the score for each complex
    for name, data in complexes.items():
        score = 0
        score_calculation = []
        for feature in data['features']:
            points = stability_scores.get(feature, 0)
            score += points
            score_calculation.append(str(points))
        
        final_scores[name] = score
        
        print(f"Analysis for {name}:")
        print(f"  - Cyclometalating Ligand: {data['cyclometalating_ligand']}")
        print(f"  - Ancillary Ligand: {data['ancillary_ligand']}")
        if not score_calculation:
            score_calculation.append("0")
        
        # Output the "equation" as requested
        print(f"  - Stability Score Calculation: {' + '.join(score_calculation)} = {score}")
        print("-" * 60)

    # 4. Determine the most stable complex and provide a conclusion
    most_stable_complex = max(final_scores, key=final_scores.get)

    print("\nConclusion:")
    print(f"Complex 3 has the highest stability score of {final_scores[most_stable_complex]}.")
    print("It effectively combines two key stabilizing strategies:")
    print("  1. Fluorination to strengthen the Ir-C bond against chemical degradation.")
    print("  2. Bulky tert-butyl groups to prevent intermolecular quenching and improve film morphology.")
    print("\nTherefore, LECs based on Complex 3 are expected to be the most stable.")

analyze_emitter_stability()
<<<C>>>