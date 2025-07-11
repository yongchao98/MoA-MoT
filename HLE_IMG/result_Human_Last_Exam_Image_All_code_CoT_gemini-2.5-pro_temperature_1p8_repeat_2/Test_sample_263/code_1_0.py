def analyze_lec_stability():
    """
    Analyzes the expected stability of LECs based on three Iridium complexes.

    This function simulates the chemical reasoning by assigning stability scores
    to different molecular features.
    """

    # Define the structural features of each complex relevant to stability
    complexes = {
        'Complex 1': {
            'description': 'Baseline complex with standard ppy and bpy ligands.',
            'features': []
        },
        'Complex 2': {
            'description': 'Complex with a large, planar ancillary ligand.',
            'features': ['large_planar_ligand']
        },
        'Complex 3': {
            'description': 'Complex with fluorinated and bulky ligands.',
            'features': ['fluorination', 'bulky_groups']
        }
    }

    # Assign scores to each feature.
    # Positive scores enhance stability, negative scores may decrease it.
    stability_scores = {
        'fluorination': 1,  # Strengthens Ir-C bonds, enhancing stability.
        'bulky_groups': 1,  # Provide steric protection, enhancing stability.
        'large_planar_ligand': -1 # Can lead to aggregation, a potential degradation pathway.
    }

    print("Analyzing the expected stability of each complex...\n")

    # Calculate and store the stability score for each complex
    ranked_complexes = {}
    for name, data in complexes.items():
        score = 0
        print(f"--- {name} ---")
        print(data['description'])
        print("Features contributing to stability:")
        if not data['features']:
            print("- None (Baseline)")
        else:
            for feature in data['features']:
                feature_score = stability_scores.get(feature, 0)
                score += feature_score
                print(f"- {feature.replace('_', ' ')} (Score: {feature_score})")
        
        ranked_complexes[name] = score
        print(f"Total Stability Score for {name}: {score}\n")

    # Determine the most stable complex
    most_stable_complex = max(ranked_complexes, key=ranked_complexes.get)
    
    print("--- Conclusion ---")
    print(f"The complex with the highest stability score is {most_stable_complex}.")
    print("Reasoning: Complex 3 possesses two key stabilizing features:")
    print("1. Fluorination of the main ligands, which strengthens the Ir-C bonds.")
    print("2. Bulky tert-butyl groups on the ancillary ligand, which provide steric protection.")
    print("These modifications are known to significantly increase the operational lifetime of light-emitting devices.")
    print("\nTherefore, LECs based on complex 3 are expected to result in more stable devices.")

# Run the analysis
analyze_lec_stability()