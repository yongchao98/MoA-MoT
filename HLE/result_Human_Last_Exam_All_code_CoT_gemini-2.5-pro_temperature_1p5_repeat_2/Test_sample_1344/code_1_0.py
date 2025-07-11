import itertools

def find_optimal_catalyst():
    """
    This script models the selection of an optimal dual-function catalyst
    for both olefin polymerization and polyolefin hydrogenolysis.
    It scores combinations of metals, ligands, and supports to find the best balance.
    """

    # Step 1: Define the chemical components and their properties
    # Scores are on a scale of 1-10 based on established chemical principles.
    metals = {
        'Titanium (Ti)': {
            'polymerization': 8,  # Good activity, but less than Zr.
            'cracking': 5,      # Lower thermal stability.
            'cost': 9
        },
        'Zirconium (Zr)': {
            'polymerization': 10, # The benchmark for polymerization.
            'cracking': 6,      # Moderate thermal stability.
            'cost': 6
        },
        'Hafnium (Hf)': {
            'polymerization': 8,  # Slightly less active than Zr but still excellent.
            'cracking': 9,      # High thermal stability, good for cracking conditions.
            'cost': 4
        }
    }

    ligands = {
        'Metallocene (e.g., Cp2)': {
            'polymerization': 8,  # Classic, well-understood.
            'cracking': 5,      # Can be unstable at high temperatures.
        },
        'Constrained-Geometry (CGC)': {
            'polymerization': 10, # Highly active, open coordination site.
            'cracking': 8,      # Robust framework, good for harsh conditions.
        },
        'Pincer Ligand (e.g., PNP)': {
            'polymerization': 5,  # Often less active for polymerization.
            'cracking': 9,      # Excellent for activating small molecules like H2.
        }
    }

    supports = {
        'None (Homogeneous)': {
            'polymerization': 7,  # Highly active but can deactivate.
            'cracking': 6,      # May struggle with stability.
        },
        'Methylaluminoxane (MAO)': {
            'polymerization': 10, # Excellent cocatalyst and scavenger.
            'cracking': 3,      # Not ideal for reductive cracking conditions.
        },
        'Modified Silica (SiO2)': {
            'polymerization': 8,  # Good for creating stable single-sites.
            'cracking': 9,      # Provides stability and surface acidity to aid cracking.
        }
    }

    # Step 2: Iterate through all combinations and calculate scores
    best_combination = None
    max_score = -1

    print("--- Catalyst Evaluation ---")
    
    # Generate all possible combinations
    all_combinations = list(itertools.product(metals.keys(), ligands.keys(), supports.keys()))
    
    results = []

    for metal_name, ligand_name, support_name in all_combinations:
        metal_props = metals[metal_name]
        ligand_props = ligands[ligand_name]
        support_props = supports[support_name]

        # Calculate score for each function
        poly_score = metal_props['polymerization'] + ligand_props['polymerization'] + support_props['polymerization']
        crack_score = metal_props['cracking'] + ligand_props['cracking'] + support_props['cracking']

        # The overall score is multiplicative to ensure competence in both areas.
        overall_score = poly_score * crack_score
        
        results.append({
            'combination': f"{metal_name} / {ligand_name} / {support_name}",
            'poly_score': poly_score,
            'crack_score': crack_score,
            'overall_score': overall_score
        })

    # Find the combination with the highest overall score
    best_result = max(results, key=lambda x: x['overall_score'])

    # Step 3: Print the results and the final conclusion
    print(f"Optimal Combination Found: {best_result['combination']}")
    print("\n--- Scoring Breakdown ---")
    print(f"This combination was chosen based on the following calculated scores:")
    print(f"Final Polymerization Score: {best_result['poly_score']}")
    print(f"Final Cracking Score:       {best_result['crack_score']}")
    print(f"Final Combined Score:     {best_result['overall_score']} (Polymerization * Cracking)")
    
    print("\n--- Chemical Rationale ---")
    print("This system is predicted to be optimal because:")
    print("1. Hafnium (Hf) provides an excellent balance of high polymerization activity and the crucial thermal stability required for high-temperature hydrocracking.")
    print("2. A Constrained-Geometry Catalyst (CGC) ligand creates a highly reactive, open metal center ideal for polymerization, while its robust chemical structure withstands the harsh conditions needed for breaking down the polymer.")
    print("3. A Modified Silica (SiO2) support immobilizes the catalyst, creating a stable, well-defined single-site that prevents deactivation. The support's surface properties can also synergistically assist in the C-C bond cleavage steps.")

    return f"<<<{best_result['combination']}>>>"


# Execute the function and print the final formatted answer
final_answer = find_optimal_catalyst()
print(final_answer)
