import itertools

def evaluate_bifunctional_catalysts():
    """
    This script models the search for an optimal single-site catalyst
    for both olefin polymerization and polyolefin hydrogenolysis.
    It ranks combinations of metals, ligands, and supports based on a
    simplified scoring system derived from known chemical principles.
    """

    # 1. Define components and their scores based on chemical suitability
    # Scores are out of 10 for each category.
    components = {
        'metal': {
            'Titanium (Ti)': {'poly_score': 7, 'hydro_score': 6},
            'Zirconium (Zr)': {'poly_score': 9, 'hydro_score': 10},
            'Hafnium (Hf)': {'poly_score': 10, 'hydro_score': 8},
        },
        'ligand': {
            'Metallocene': {'poly_score': 10, 'hydro_score': 4},
            'Constrained-Geometry (CGC)': {'poly_score': 8, 'hydro_score': 10},
            'Phenoxy-imine (FI)': {'poly_score': 9, 'hydro_score': 6},
        },
        'support': {
            'None (Homogeneous)': {'poly_score': 10, 'hydro_score': 3},
            'Silica (SiO2)': {'poly_score': 7, 'hydro_score': 9},
            'Zirconia (ZrO2)': {'poly_score': 8, 'hydro_score': 10},
        }
    }

    # 2. Iterate through all combinations and calculate scores
    best_combination = None
    max_bifunctional_score = -1

    metal_options = components['metal'].keys()
    ligand_options = components['ligand'].keys()
    support_options = components['support'].keys()

    for metal, ligand, support in itertools.product(metal_options, ligand_options, support_options):
        metal_scores = components['metal'][metal]
        ligand_scores = components['ligand'][ligand]
        support_scores = components['support'][support]

        # Calculate scores for each function
        total_poly_score = metal_scores['poly_score'] + ligand_scores['poly_score'] + support_scores['poly_score']
        total_hydro_score = metal_scores['hydro_score'] + ligand_scores['hydro_score'] + support_scores['hydro_score']
        
        # The bifunctional score is the sum of the two individual scores.
        # This rewards combinations that are strong in both categories.
        # The equation is: Bifunctional Score = Polymerization Score + Hydrogenolysis Score
        bifunctional_score = total_poly_score + total_hydro_score

        if bifunctional_score > max_bifunctional_score:
            max_bifunctional_score = bifunctional_score
            best_combination = {
                'Metal': metal,
                'Ligand': ligand,
                'Support': support,
                'Polymerization Score': total_poly_score,
                'Hydrogenolysis Score': total_hydro_score,
                'Bifunctional Score': bifunctional_score
            }

    # 3. Print the results for the optimal combination
    print("--- Optimal Bifunctional Catalyst Analysis ---")
    if best_combination:
        print(f"\nThe optimal combination found by the model is:")
        print(f"  - Metal:   {best_combination['Metal']}")
        print(f"  - Ligand:  {best_combination['Ligand']}")
        print(f"  - Support: {best_combination['Support']}")
        
        print("\n--- Scoring Breakdown ---")
        print(f"Suitability for Polymerization: {best_combination['Polymerization Score']} / 30")
        print(f"Suitability for Hydrogenolysis: {best_combination['Hydrogenolysis Score']} / 30")
        
        print("\nFinal Scoring Equation and Result:")
        # The prompt requires printing each number in the final equation.
        poly_score_val = best_combination['Polymerization Score']
        hydro_score_val = best_combination['Hydrogenolysis Score']
        final_score_val = best_combination['Bifunctional Score']
        
        print(f"Final Bifunctional Score = Polymerization Score + Hydrogenolysis Score")
        print(f"Final Score -> {final_score_val} = {poly_score_val} + {hydro_score_val}")
    else:
        print("Could not determine an optimal combination.")

# Execute the analysis
evaluate_bifunctional_catalysts()