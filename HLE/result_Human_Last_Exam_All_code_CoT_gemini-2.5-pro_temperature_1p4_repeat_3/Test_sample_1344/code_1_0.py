import itertools

def find_optimal_catalyst():
    """
    Models the search for an optimal dual-function catalyst for
    olefin polymerization and polyolefin hydrogenolysis.
    """
    # 1. Define Components and their hypothetical performance scores
    # Scores are based on general chemical principles for each task.
    metals = [
        {'name': 'Titanium (Ti)', 'poly_score': 7, 'hydro_score': 6},
        {'name': 'Zirconium (Zr)', 'poly_score': 9, 'hydro_score': 8},
        {'name': 'Hafnium (Hf)', 'poly_score': 10, 'hydro_score': 5}
    ]

    ligands = [
        {'name': 'Metallocene (e.g., Cp*_2)', 'poly_score': 9, 'hydro_score': 4},
        {'name': 'Constrained-Geometry (CGC)', 'poly_score': 8, 'hydro_score': 9},
        {'name': 'Phenoxy-imine (FI)', 'poly_score': 7, 'hydro_score': 7}
    ]

    supports = [
        {'name': 'Unsupported (Homogeneous)', 'poly_score': 0, 'hydro_score': 5},
        {'name': 'Silica (SiO2)', 'poly_score': 2, 'hydro_score': 1},
        {'name': 'MAO-on-Silica', 'poly_score': 3, 'hydro_score': 2}
    ]

    # 2. Initialize variables to track the best combination
    best_combination = None
    max_optimality_score = -1
    best_scores = {}

    # 3. Iterate through all possible combinations
    for metal, ligand, support in itertools.product(metals, ligands, supports):
        # 4. Calculate total scores for the current combination
        total_poly_score = metal['poly_score'] + ligand['poly_score'] + support['poly_score']
        total_hydro_score = metal['hydro_score'] + ligand['hydro_score'] + support['hydro_score']

        # The optimality metric is the product, rewarding balanced performance
        # A catalyst must be effective at both tasks.
        optimality_score = total_poly_score * total_hydro_score

        # 5. Check if this combination is the new best
        if optimality_score > max_optimality_score:
            max_optimality_score = optimality_score
            best_combination = {
                'Metal': metal,
                'Ligand': ligand,
                'Support': support
            }
            best_scores = {
                'poly': total_poly_score,
                'hydro': total_hydro_score
            }

    # 6. Print the results
    print("--- Catalyst Optimization Model ---")
    print("Searching for an optimal catalyst with dual functionality for polymerization and hydrogenolysis.\n")

    if best_combination:
        print("Optimal combination found:")
        print(f"  - Metal:   {best_combination['Metal']['name']}")
        print(f"  - Ligand:  {best_combination['Ligand']['name']}")
        print(f"  - Support: {best_combination['Support']['name']}\n")

        print("--- Performance Calculation ---")
        # The final "equation" showing how the scores are calculated
        print("Final Polymerization Score = "
              f"{best_combination['Metal']['poly_score']} (Metal) + "
              f"{best_combination['Ligand']['poly_score']} (Ligand) + "
              f"{best_combination['Support']['poly_score']} (Support) = {best_scores['poly']}")

        print("Final Hydrogenolysis Score = "
              f"{best_combination['Metal']['hydro_score']} (Metal) + "
              f"{best_combination['Ligand']['hydro_score']} (Ligand) + "
              f"{best_combination['Support']['hydro_score']} (Support) = {best_scores['hydro']}")

        print("\nOverall Optimality Metric (Poly Score * Hydro Score):")
        print(f"{best_scores['poly']} * {best_scores['hydro']} = {max_optimality_score}")

    else:
        print("No suitable catalyst combination could be determined.")

if __name__ == '__main__':
    find_optimal_catalyst()