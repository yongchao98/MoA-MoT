import pandas as pd

def evaluate_catalyst_systems():
    """
    Evaluates and ranks catalyst systems for both olefin polymerization and
    polyolefin hydrogenolysis based on a simplified scoring model.
    """
    # Define candidate catalyst systems
    candidates = [
        {'id': 'A', 'metal': 'Zr', 'ligand': 'ansa-Indenyl', 'support': 'None'},
        {'id': 'B', 'metal': 'Hf', 'ligand': 'CGC', 'support': 'None'},
        {'id': 'C', 'metal': 'Zr', 'ligand': 'Pincer', 'support': 'Al2O3'},
        {'id': 'D', 'metal': 'Hf', 'ligand': 'Cp*', 'support': 'ZrO2'},
        {'id': 'E', 'metal': 'Hf', 'ligand': 'CGC', 'support': 'ZrO2'},
    ]

    # Scoring matrices based on chemical principles
    # Hf/Zr are superior to Ti for polymerization; CGC ligands are state-of-the-art.
    # For hydrogenolysis (C-C cleavage), a more electrophilic center (Zr) and
    # an active support are beneficial. Pincer ligands are known for bond activation.
    # For stability, Hf is generally more robust, and bulky/chelating ligands/supports help.
    scores = {
        'polymerization': {'metal': {'Hf': 3, 'Zr': 2, 'Ti': 1},
                           'ligand': {'CGC': 3, 'ansa-Indenyl': 2, 'Cp*': 1, 'Pincer': 0}},
        'hydrogenolysis': {'metal': {'Zr': 3, 'Hf': 2, 'Ti': 1},
                           'ligand': {'Pincer': 3, 'CGC': 2, 'Cp*': 1, 'ansa-Indenyl': 1},
                           'support': {'Al2O3': 2, 'ZrO2': 2, 'None': 0}},
        'stability': {'metal': {'Hf': 3, 'Zr': 2, 'Ti': 1},
                      'ligand': {'Pincer': 3, 'Cp*': 2, 'CGC': 2, 'ansa-Indenyl': 1},
                      'support': {'Al2O3': 1, 'ZrO2': 1, 'None': 0}}
    }

    # Weights to prioritize the more difficult hydrogenolysis task
    weights = {
        'polymerization': 1.0,
        'hydrogenolysis': 1.5,
        'stability': 1.0
    }

    results = []
    for cat in candidates:
        poly_score = scores['polymerization']['metal'].get(cat['metal'], 0) + \
                     scores['polymerization']['ligand'].get(cat['ligand'], 0)

        hydro_score = scores['hydrogenolysis']['metal'].get(cat['metal'], 0) + \
                      scores['hydrogenolysis']['ligand'].get(cat['ligand'], 0) + \
                      scores['hydrogenolysis']['support'].get(cat['support'], 0)

        stab_score = scores['stability']['metal'].get(cat['metal'], 0) + \
                     scores['stability']['ligand'].get(cat['ligand'], 0) + \
                     scores['stability']['support'].get(cat['support'], 0)
        
        total_score = (poly_score * weights['polymerization'] +
                       hydro_score * weights['hydrogenolysis'] +
                       stab_score * weights['stability'])

        results.append({
            'System': f"{cat['metal']}/{cat['ligand']}/{cat['support']}",
            'Polymerization Score': poly_score,
            'Hydrogenolysis Score': hydro_score,
            'Stability Score': stab_score,
            'Final Weighted Score': round(total_score, 2)
        })

    # Determine the best catalyst
    results_df = pd.DataFrame(results)
    best_catalyst = results_df.loc[results_df['Final Weighted Score'].idxmax()]

    print("--- Catalyst Evaluation Model ---")
    print("Scoring candidates for dual-function catalysis:")
    print(results_df.to_string(index=False))
    print("\n--- Optimal Combination Identified ---")
    print(f"The best-performing system based on this model is: {best_catalyst['System']}")
    print("\nThis system provides the best balance of high activity for both reactions and overall stability.")
    print("The final scoring equation for the optimal catalyst is:")
    
    # Fulfilling the "output each number in the final equation" requirement
    poly = best_catalyst['Polymerization Score']
    hydro = best_catalyst['Hydrogenolysis Score']
    stab = best_catalyst['Stability Score']
    w_poly = weights['polymerization']
    w_hydro = weights['hydrogenolysis']
    w_stab = weights['stability']
    total = best_catalyst['Final Weighted Score']
    
    print(f"({poly} * {w_poly})_Polymerization + ({hydro} * {w_hydro})_Hydrogenolysis + ({stab} * {w_stab})_Stability = {total} Total Score")

evaluate_catalyst_systems()