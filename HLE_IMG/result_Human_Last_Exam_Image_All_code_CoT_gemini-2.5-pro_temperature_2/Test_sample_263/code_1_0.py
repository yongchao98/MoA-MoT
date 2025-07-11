def analyze_lec_stability():
    """
    Analyzes the expected stability of LECs based on three Ir(III) complexes.
    The analysis is based on well-known molecular design principles for improving
    the stability of organic electronic materials.
    """

    # Define the complexes and their key structural features
    complexes = {
        'Complex 1': {
            'cyclometalating_ligand': 'ppy (phenylpyridine)',
            'diimine_ligand': 'bpy (bipyridine)',
            'features': []
        },
        'Complex 2': {
            'cyclometalating_ligand': 'ppy (phenylpyridine)',
            'diimine_ligand': 'phen-imi (phenanthroline-imidazole derivative)',
            'features': ['extended_pi_system_on_diimine']
        },
        'Complex 3': {
            'cyclometalating_ligand': 'dfppy (difluoro-phenylpyridine)',
            'diimine_ligand': 'dtbbpy (di-tert-butyl-bipyridine)',
            'features': ['fluorination', 'bulky_substituents']
        }
    }

    # Assign stability scores based on features
    # Base score = 1
    # Fluorination = +1 (enhances electrochemical stability against oxidation)
    # Bulky substituents = +1 (improves morphological stability and steric protection)
    # Extended pi system = 0 (effect on stability is not clearly predictable as positive)
    
    stability_scores = {}
    print("Evaluating the stability of each complex based on its structure:\n")

    for name, data in complexes.items():
        score = 1  # Base score
        print(f"--- Analyzing {name} ---")
        print(f"  Ligands: {data['cyclometalating_ligand']} and {data['diimine_ligand']}")
        
        reasoning = []
        if 'fluorination' in data['features']:
            score += 1
            reasoning.append("  - Contains fluorinated ligands, which increases resistance to electrochemical oxidation (+1 score).")
        if 'bulky_substituents' in data['features']:
            score += 1
            reasoning.append("  - Contains bulky tert-butyl groups, which provide steric protection and improve morphological stability (+1 score).")
        
        if not reasoning:
            reasoning.append("  - This is the reference complex with no special stability-enhancing modifications (base score).")

        for reason in reasoning:
            print(reason)
            
        stability_scores[name] = score
        print(f"  Final Stability Score for {name}: {score}\n")

    # Determine the most stable complex
    most_stable_complex = max(stability_scores, key=stability_scores.get)
    max_score = stability_scores[most_stable_complex]

    print("--- Conclusion ---")
    print("Comparison of final stability scores:")
    for name, score in stability_scores.items():
        print(f"{name}: {score}")

    print(f"\n{most_stable_complex} has the highest stability score of {max_score}.")
    print("Therefore, LECs based on Complex 3 are expected to be the most stable due to the combined positive effects of ligand fluorination and bulky substituents.")

# Execute the analysis
analyze_lec_stability()