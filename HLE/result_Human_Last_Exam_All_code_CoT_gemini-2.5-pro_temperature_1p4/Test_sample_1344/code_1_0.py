import itertools

def find_optimal_catalyst():
    """
    A simplified model to find an optimal dual-function catalyst for
    polyolefin polymerization and hydrogenolysis.
    """
    # 1. Define Candidate Components with Hypothetical Scores
    # Scores are for 'polymerization_activity' and 'breakdown_activity'
    metals = {
        "Ti": {"poly": 6, "breakdown": 9},  # Titanium: Potentially better for breakdown
        "Zr": {"poly": 9, "breakdown": 6},  # Zirconium: A classic for polymerization
        "Hf": {"poly": 8, "breakdown": 5},  # Hafnium: Similar to Zirconium
    }

    ligands = {
        "Metallocene": {"poly": 9, "breakdown": 3}, # Great for polymerization, but sterically hindered for breakdown
        "Constrained-Geometry": {"poly": 7, "breakdown": 8}, # Open site is good for both functions
        "Phenoxy-imine": {"poly": 8, "breakdown": 6}, # A versatile modern ligand
    }

    supports = {
        "None": {"factor": 1.0}, # No support, baseline performance
        "Silica (SiO2)": {"factor": 1.1}, # Inert support, slight stability boost
        "MAO (activator)": {"factor": 1.3}, # Co-catalyst, enhances overall activity
    }

    best_catalyst = None
    max_balanced_score = -1
    best_scores = {}

    # 2. Iterate through all combinations
    catalyst_components = [metals.keys(), ligands.keys(), supports.keys()]
    for combo in itertools.product(*catalyst_components):
        metal_name, ligand_name, support_name = combo

        # Get the base scores from the dictionaries
        metal_scores = metals[metal_name]
        ligand_scores = ligands[ligand_name]
        support_factor = supports[support_name]['factor']

        # 3. Calculate scores for each function
        poly_score = (metal_scores['poly'] + ligand_scores['poly']) * support_factor
        breakdown_score = (metal_scores['breakdown'] + ligand_scores['breakdown']) * support_factor

        # 4. Calculate the balanced performance score (product of the two activities)
        balanced_score = poly_score * breakdown_score

        # 5. Check if this is the best catalyst so far
        if balanced_score > max_balanced_score:
            max_balanced_score = balanced_score
            best_catalyst = combo
            best_scores = {
                'poly': poly_score,
                'breakdown': breakdown_score,
                'balanced': balanced_score,
                'details': {
                    'metal': metal_scores,
                    'ligand': ligand_scores,
                    'support_factor': support_factor
                }
            }

    # Print the results for the best combination found
    print("--- Catalyst Optimization Model ---")
    if best_catalyst:
        metal, ligand, support = best_catalyst
        details = best_scores['details']

        print(f"\nOptimal Combination Found: {metal} center, {ligand} ligand, with {support} support/activator.\n")

        # Fulfilling the request to show the numbers in the final equation
        print("Performance Calculation:")

        poly_m_score = details['metal']['poly']
        poly_l_score = details['ligand']['poly']
        supp_factor = details['support_factor']
        print(f"  Polymerization Score: ({poly_m_score} [metal] + {poly_l_score} [ligand]) * {supp_factor} [support] = {best_scores['poly']:.2f}")

        break_m_score = details['metal']['breakdown']
        break_l_score = details['ligand']['breakdown']
        print(f"  Breakdown Score:      ({break_m_score} [metal] + {break_l_score} [ligand]) * {supp_factor} [support] = {best_scores['breakdown']:.2f}")

        print(f"\n  Final Balanced Score (Poly * Breakdown): {best_scores['poly']:.2f} * {best_scores['breakdown']:.2f} = {best_scores['balanced']:.2f}")
    else:
        print("Could not determine an optimal catalyst.")

# Run the simulation
find_optimal_catalyst()