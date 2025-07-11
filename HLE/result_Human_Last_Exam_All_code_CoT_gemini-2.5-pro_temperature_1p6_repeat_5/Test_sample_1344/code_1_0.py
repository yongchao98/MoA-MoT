import itertools

def find_optimal_catalyst():
    """
    This script models the selection of an optimal bifunctional catalyst
    for both olefin polymerization and polyolefin hydrogenolysis.
    It does so by evaluating combinations of metals, ligands, and supports
    based on a simplified scoring system.
    """
    # 1. Define Catalyst Components and their hypothetical scores
    # Scores are on a scale of 1-10 for inherent capability.
    metals = {
        'Titanium (Ti)': {'poly_score': 7, 'hydro_score': 6},
        'Zirconium (Zr)': {'poly_score': 9, 'hydro_score': 8},
        'Hafnium (Hf)': {'poly_score': 10, 'hydro_score': 7},
    }

    # Ligands influence steric and electronic properties.
    ligands = {
        'Cp* (Pentamethylcyclopentadienyl)': {'poly_score': 8, 'hydro_score': 5},
        'CGC (Constrained-Geometry)': {'poly_score': 10, 'hydro_score': 6},
        'PNP Pincer': {'poly_score': 4, 'hydro_score': 9},
    }

    # Supports can immobilize the catalyst and provide tandem functionality.
    # Scores here are bonuses added to the main scores.
    supports = {
        'None (Homogeneous)': {'poly_bonus': 0, 'hydro_bonus': 0},
        'Silica (SiO2)': {'poly_bonus': 1, 'hydro_bonus': 2},
        'Alumina (Al2O3)': {'poly_bonus': 0, 'hydro_bonus': 3},
        'MAO (as support/activator)': {'poly_bonus': 2, 'hydro_bonus': 0},
    }

    best_combination = None
    max_bifunctional_score = -1

    print("Evaluating potential catalyst combinations...\n")

    # 2. Iterate through all combinations
    for metal_name, metal_scores in metals.items():
        for ligand_name, ligand_scores in ligands.items():
            for support_name, support_scores in supports.items():
                
                # 3. Calculate total scores for each function
                total_poly_score = metal_scores['poly_score'] + ligand_scores['poly_score'] + support_scores['poly_bonus']
                total_hydro_score = metal_scores['hydro_score'] + ligand_scores['hydro_score'] + support_scores['hydro_bonus']
                
                # 4. Calculate the 'bifunctionality score'
                # This metric rewards catalysts that are effective at both reactions.
                bifunctional_score = total_poly_score * total_hydro_score

                # 5. Check if this is the best combination so far
                if bifunctional_score > max_bifunctional_score:
                    max_bifunctional_score = bifunctional_score
                    best_combination = {
                        'Metal': metal_name,
                        'Ligand': ligand_name,
                        'Support': support_name,
                        'Polymerization Score': total_poly_score,
                        'Hydrogenolysis Score': total_hydro_score,
                        'Bifunctionality Score': bifunctional_score
                    }
    
    # 6. Output the results
    if best_combination:
        print("--- Optimal Catalyst Combination Found ---\n")
        print(f"Metal: {best_combination['Metal']}")
        print(f"Ligand: {best_combination['Ligand']}")
        print(f"Support: {best_combination['Support']}\n")
        print(f"This combination was chosen because it provides the best-balanced performance for both required tasks based on the model.")
        print("The final equation for its score is:")
        print(f"Total Polymerization Score = {best_combination['Polymerization Score']}")
        print(f"Total Hydrogenolysis Score = {best_combination['Hydrogenolysis Score']}")
        print(f"Final Bifunctionality Score (Poly Score * Hydro Score) = {best_combination['Bifunctionality Score']}\n")

        # Construct the final answer string
        final_answer = (
            f"The optimal combination is {best_combination['Metal']} "
            f"with a {best_combination['Ligand']} ligand, "
            f"supported on {best_combination['Support']}."
        )
        print(f"<<<{final_answer}>>>")
    else:
        print("Could not determine an optimal catalyst.")

if __name__ == '__main__':
    find_optimal_catalyst()