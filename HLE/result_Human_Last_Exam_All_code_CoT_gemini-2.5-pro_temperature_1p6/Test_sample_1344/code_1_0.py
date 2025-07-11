import itertools

def find_optimal_catalyst():
    """
    This function models the search for an optimal dual-function catalyst.
    It uses hypothetical scores to evaluate combinations of metals, ligands, and supports
    for both polymerization and hydrogenolysis capabilities.
    """

    # 1. Define components and their hypothetical performance scores (out of 10)
    # These scores are based on general chemical principles for illustration.
    # A higher score indicates better performance for that specific task.
    components = {
        'metal': {
            'Ti': {'poly': 7, 'hydro': 5},
            'Zr': {'poly': 9, 'hydro': 8},
            'Hf': {'poly': 8, 'hydro': 9},
        },
        'ligand': {
            # Metallocene (Cp*2): Excellent for polymerization
            'Cp*2': {'poly': 8, 'hydro': 3},
            # Constrained-Geometry Catalyst (CGC): Balances polymerization and stability
            'CGC': {'poly': 9, 'hydro': 5},
            # Pincer Ligand (PNP): Often used in hydrogenation/dehydrogenation reactions
            'PNP': {'poly': 4, 'hydro': 9},
        },
        'support': {
            # Homogeneous (no support)
            'None': {'poly': 0, 'hydro': 0},
            # Common inert support
            'SiO2': {'poly': 0, 'hydro': 1},
            # Support with Lewis acidic properties
            'Al2O3': {'poly': 1, 'hydro': 2},
            # Tunable, porous support that can aid in tandem catalysis
            'MOF': {'poly': 3, 'hydro': 4},
        }
    }

    # 2. Define the weights for each catalytic task. We want a balanced catalyst.
    poly_weight = 0.5
    hydro_weight = 0.5

    # 3. Initialize variables to track the best combination found
    best_combination = None
    max_score = -1
    best_scores_detail = {}

    # 4. Generate all possible combinations
    metals = components['metal'].keys()
    ligands = components['ligand'].keys()
    supports = components['support'].keys()
    all_combinations = list(itertools.product(metals, ligands, supports))

    # 5. Iterate through each combination and calculate its score
    for combo in all_combinations:
        metal, ligand, support = combo

        # Calculate the total score for each task by summing component scores
        poly_score = (components['metal'][metal]['poly'] +
                      components['ligand'][ligand]['poly'] +
                      components['support'][support]['poly'])

        hydro_score = (components['metal'][metal]['hydro'] +
                       components['ligand'][ligand]['hydro'] +
                       components['support'][support]['hydro'])

        # Calculate the final weighted score
        total_score = (poly_score * poly_weight) + (hydro_score * hydro_weight)

        # If this combination is the best so far, save it
        if total_score > max_score:
            max_score = total_score
            best_combination = combo
            best_scores_detail = {
                'poly': poly_score,
                'hydro': hydro_score,
                'total': round(total_score, 2)
            }

    # 6. Print the results
    if best_combination:
        print("--- Catalyst Optimization Model ---")
        print("Based on the model, the optimal combination is:")
        print(f"  - Group IV Metal: {best_combination[0]}")
        print(f"  - Ligand: {best_combination[1]}")
        print(f"  - Support: {best_combination[2]}")
        print("\nThis combination provides the best balance of modeled capabilities:")
        print(f"  - Polymerization Score: {best_scores_detail['poly']}")
        print(f"  - Hydrogenolysis Score: {best_scores_detail['hydro']}")

        print("\nFinal Score Calculation:")
        final_poly_score = best_scores_detail['poly']
        final_hydro_score = best_scores_detail['hydro']
        final_total_score = best_scores_detail['total']
        
        # Output the numbers used in the final equation
        print(f"Total Score = (Polymerization Score * Weight) + (Hydrogenolysis Score * Weight)")
        print(f"{final_total_score} = ({final_poly_score} * {poly_weight}) + ({final_hydro_score} * {hydro_weight})")

    else:
        print("Could not determine an optimal combination.")

# Run the analysis
find_optimal_catalyst()
<<<Hf>>>