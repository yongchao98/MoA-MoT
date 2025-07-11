import pandas as pd

def analyze_lec_stability():
    """
    Analyzes the expected stability of three Ir(III) complexes in LECs
    based on their structural features.
    """
    
    # Define stability scoring rules based on common molecular design principles.
    # A higher score indicates a higher expected operational stability.
    scoring_rules = {
        'base_score': 5,  # Baseline for a heteroleptic Ir(III) complex
        'fluorination_per_ligand': 2, # Fluorination enhances electrochemical stability
        'bulky_groups': 3, # Bulky groups prevent aggregation, improving morphological stability
        'large_planar_ligand': -2, # Large planar ligands can promote aggregation, reducing stability
        'unsubstituted_ligand': 0 # No significant impact on stability relative to baseline
    }

    # Define the features for each complex based on the provided image
    complexes = {
        'Complex 1': {
            'C^N_ligand': 'unsubstituted', 'N^N_ligand': 'unsubstituted',
            'description': "[Ir(ppy)2(bpy)]+"
        },
        'Complex 2': {
            'C^N_ligand': 'unsubstituted', 'N^N_ligand': 'large_planar',
            'description': "[Ir(ppy)2(large_planar_N^N)]+"
        },
        'Complex 3': {
            'C^N_ligand': 'fluorinated', 'N^N_ligand': 'bulky_groups',
            'description': "[Ir(dfppy)2(dtbpy)]+"
        }
    }
    
    results = []

    print("Analyzing LEC Stability based on Molecular Structure\n")
    print("Stability scoring system:")
    for rule, points in scoring_rules.items():
        print(f"- {rule.replace('_', ' ').capitalize()}: {points} points")
    print("-" * 50)
    
    # Calculate stability score for each complex
    for name, features in complexes.items():
        score = scoring_rules['base_score']
        calculation_steps = [f"{scoring_rules['base_score']} (base score)"]
        
        # Score C^N ligands (there are two)
        if features['C^N_ligand'] == 'fluorinated':
            score += 2 * scoring_rules['fluorination_per_ligand']
            calculation_steps.append(f"2 * {scoring_rules['fluorination_per_ligand']} (fluorinated C^N ligands)")
        
        # Score N^N ligand
        if features['N^N_ligand'] == 'bulky_groups':
            score += scoring_rules['bulky_groups']
            calculation_steps.append(f"{scoring_rules['bulky_groups']} (bulky N^N ligand)")
        elif features['N^N_ligand'] == 'large_planar':
            score += scoring_rules['large_planar_ligand']
            calculation_steps.append(f"{scoring_rules['large_planar_ligand']} (large planar N^N ligand)")
        
        results.append({'Complex': name, 'Score': score, 'Calculation': calculation_steps})

    # Print results and reasoning
    max_score = -1
    most_stable_complex = None
    
    for res in results:
        print(f"\nAnalysis for {res['Complex']} ({features['description']}):")
        
        # Print the reasoning from the scoring
        if res['Complex'] == 'Complex 1':
            print("Reasoning: This is a standard benchmark complex with no special stabilizing groups.")
        elif res['Complex'] == 'Complex 2':
            print("Reasoning: This complex has a large, planar N^N ligand, which is likely to cause aggregation and reduce device stability.")
        elif res['Complex'] == 'Complex 3':
            print("Reasoning: This complex has two key stabilizing features:")
            print("1. Fluorination of the C^N ligands increases electrochemical stability.")
            print("2. Bulky tert-butyl groups on the N^N ligand prevent aggregation and improve morphological stability.")
        
        # Print the score calculation
        calc_str = " + ".join(calculation_steps).replace('+ -', '- ')
        print(f"Stability Score Calculation: {res['Score']} = " + " + ".join(res['Calculation']).replace("+ -", "- "))
        
        if res['Score'] > max_score:
            max_score = res['Score']
            most_stable_complex = res['Complex']

    print("\n" + "-" * 50)
    print(f"Conclusion: {most_stable_complex} has the highest stability score ({max_score}).")
    print("Therefore, LECs based on this complex are expected to be the most stable.")

analyze_lec_stability()