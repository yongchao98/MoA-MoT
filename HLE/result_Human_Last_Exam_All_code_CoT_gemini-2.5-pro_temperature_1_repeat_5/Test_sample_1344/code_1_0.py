import itertools

def find_optimal_catalyst():
    """
    This function searches for an optimal single-site catalyst combination
    for both olefin polymerization and polyolefin hydrogenolysis.
    It uses a simplified knowledge base and scoring system to evaluate candidates.
    """

    # 1. Knowledge Base: Define components and their scores based on chemical literature.
    # Scores are on a scale of 1-10.
    metals = {
        'Zr (Zirconium)': {'polymerization': 9, 'hydrogenolysis': 9},
        'Hf (Hafnium)':   {'polymerization': 8, 'hydrogenolysis': 7},
        'Ti (Titanium)':  {'polymerization': 7, 'hydrogenolysis': 6}
    }

    ligands = {
        'Constrained-Geometry (e.g., CGC)': {
            'description': 'Excellent for incorporating comonomers, high activity.',
            'polymerization': 10,
            'hydrogenolysis': 7
        },
        'Metallocene (e.g., Cp*2)': {
            'description': 'Classic, highly active polymerization catalyst.',
            'polymerization': 9,
            'hydrogenolysis': 6
        },
        'Non-metallocene Pincer (e.g., PNP)': {
            'description': 'Highly tunable for various catalytic cycles.',
            'polymerization': 6,
            'hydrogenolysis': 8
        }
    }

    supports = {
        'Sulfonated Zirconia (ZrO2/SO4)': {
            'description': 'Creates strong acid sites, ideal for C-C cleavage.',
            'polymerization': 7, # Provides good site isolation
            'hydrogenolysis': 10 # Excellent for activating C-C bonds
        },
        'Silica (SiO2)': {
            'description': 'Common, inert support for heterogenization.',
            'polymerization': 8,
            'hydrogenolysis': 4
        },
        'Unsupported (Homogeneous)': {
            'description': 'High intrinsic activity but impractical for flow processes.',
            'polymerization': 6, # Prone to deactivation
            'hydrogenolysis': 5
        }
    }

    # 2. Search for the best combination
    best_combination = None
    highest_score = -1

    for metal_name, metal_scores in metals.items():
        for ligand_name, ligand_data in ligands.items():
            for support_name, support_scores in supports.items():
                
                # Calculate individual activity scores for the combination
                total_polymerization_score = metal_scores['polymerization'] + ligand_data['polymerization'] + support_scores['polymerization']
                total_hydrogenolysis_score = metal_scores['hydrogenolysis'] + ligand_data['hydrogenolysis'] + support_scores['hydrogenolysis']

                # The Dual-Function Score is the product of the two activities.
                # This ensures the catalyst is strong in BOTH areas.
                dual_function_score = total_polymerization_score * total_hydrogenolysis_score
                
                if dual_function_score > highest_score:
                    highest_score = dual_function_score
                    best_combination = {
                        'Metal': metal_name,
                        'Ligand': ligand_name,
                        'Support': support_name,
                        'Total Polymerization Score': total_polymerization_score,
                        'Total Hydrogenolysis Score': total_hydrogenolysis_score,
                        'Final Dual-Function Score': dual_function_score,
                        'Ligand Description': ligand_data['description'],
                        'Support Description': support_scores['description']
                    }

    # 3. Output the result
    print("--- Catalyst Optimization Search Results ---")
    if best_combination:
        print(f"Optimal Combination Found:\n")
        print(f"  Metal:   {best_combination['Metal']}")
        print(f"  Ligand:  {best_combination['Ligand']}")
        print(f"  Support: {best_combination['Support']}")
        print("\n--- Justification & Scores ---")
        print(f"This combination provides the best balance between the two opposing functions.")
        print(f"Polymerization Score Breakdown:")
        print(f"  Metal ({metals[best_combination['Metal']]['polymerization']}) + Ligand ({ligands[best_combination['Ligand']]['polymerization']}) + Support ({supports[best_combination['Support']]['polymerization']}) = {best_combination['Total Polymerization Score']}")
        print(f"Hydrogenolysis Score Breakdown:")
        print(f"  Metal ({metals[best_combination['Metal']]['hydrogenolysis']}) + Ligand ({ligands[best_combination['Ligand']]['hydrogenolysis']}) + Support ({supports[best_combination['Support']]['hydrogenolysis']}) = {best_combination['Total Hydrogenolysis Score']}")
        print("\n--- Final Calculation ---")
        print("Final Equation: (Total Polymerization Score) * (Total Hydrogenolysis Score)")
        print(f"Final Result: {best_combination['Total Polymerization Score']} * {best_combination['Total Hydrogenolysis Score']} = {best_combination['Final Dual-Function Score']}\n")
        
        # Format the final answer string as requested
        final_answer_str = f"Metal: {best_combination['Metal']}, Ligand: {best_combination['Ligand']}, Support: {best_combination['Support']}"
        print(f"<<<{final_answer_str}>>>")

    else:
        print("Could not determine an optimal combination.")

# Run the simulation
find_optimal_catalyst()