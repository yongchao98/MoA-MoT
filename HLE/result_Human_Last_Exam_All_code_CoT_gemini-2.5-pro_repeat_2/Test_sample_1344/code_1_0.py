import textwrap

def find_optimal_catalyst():
    """
    Simulates a search for an optimal bifunctional single-site catalyst
    for both olefin polymerization and polyolefin hydrogenolysis.
    """

    # Define candidate components based on established catalysis research.
    # Scores are conceptual ratings (1-10) for each task.
    metals = {
        'Zirconium (Zr)': {'polymerization': 9, 'hydrogenolysis': 7},
        'Hafnium (Hf)': {'polymerization': 8, 'hydrogenolysis': 6},
        'Titanium (Ti)': {'polymerization': 7, 'hydrogenolysis': 5}
    }

    ligands = {
        'Pyridyl-amine': {
            'polymerization': 9,
            'hydrogenolysis': 8,
            'notes': "Post-metallocene ligand; highly tunable electronically and sterically, promoting high activity."
        },
        'Constrained-Geometry (CGC)': {
            'polymerization': 10,
            'hydrogenolysis': 6,
            'notes': "Excellent for polymerization, but less inherently suited for C-C cleavage compared to more tunable systems."
        },
        'Bis(cyclopentadienyl) [Cp2]': {
            'polymerization': 7,
            'hydrogenolysis': 3,
            'notes': "Classic metallocene; less active and stable for the demanding degradation reaction."
        }
    }

    supports = {
        'Sulfated Zirconia (S-ZrO2)': {
            'synergy_score': 9,
            'notes': "A solid acid support that actively assists in C-C bond cleavage (cracking), creating a powerful bifunctional system with the metal center."
        },
        'Inert Silica (SiO2)': {
            'synergy_score': 5,
            'notes': "Provides site isolation and improves handling, but offers no chemical assistance for degradation."
        },
        'None (Homogeneous)': {
            'synergy_score': 2,
            'notes': "Lacks the synergistic effect of an active support and complicates catalyst recovery."
        }
    }

    # Evaluate all combinations to find the best one
    best_combination = None
    max_score = -1
    best_rationale = []

    for metal_name, metal_scores in metals.items():
        for ligand_name, ligand_scores in ligands.items():
            for support_name, support_scores in supports.items():
                # Calculate a score based on performance in both reactions and support synergy
                # A balanced score is preferred, so we multiply the components.
                score = (metal_scores['polymerization'] + ligand_scores['polymerization']) * \
                        (metal_scores['hydrogenolysis'] + ligand_scores['hydrogenolysis']) * \
                        support_scores['synergy_score']

                if score > max_score:
                    max_score = score
                    best_combination = {
                        "Metal": metal_name,
                        "Ligand": ligand_name,
                        "Support": support_name
                    }
                    best_rationale = [
                        f"The {metal_name} center provides an excellent balance of high activity for polymerization and reactivity towards hydrogen for the degradation step.",
                        f"The {ligand_name} ligand is a highly tunable post-metallocene system. Its electronic properties can be adjusted to create a catalyst that is both an extremely fast polymerization agent and is capable of activating the polymer C-C bonds for hydrogenolysis.",
                        f"The {support_name} support is critical. As a solid acid, it works synergistically with the single-site metal center. The metal site handles the interaction with hydrogen and the polymer chain end, while the acidic nature of the support facilitates the cleavage of the strong C-C bonds in the polymer backbone."
                    ]

    # Print the final results
    print("--- Optimal Bifunctional Catalyst Analysis ---")
    if best_combination:
        print("\nAn optimal combination for both polymerization and degradation is proposed:")
        print(f"\n  Metal:   {best_combination['Metal']}")
        print(f"  Ligand:  {best_combination['Ligand']}")
        print(f"  Support: {best_combination['Support']}")
        print("\nRationale for this combination:")
        for line in best_rationale:
            # Wrap text for better readability
            wrapped_text = textwrap.fill(line, width=80, initial_indent="- ", subsequent_indent="  ")
            print(wrapped_text)
    else:
        print("Could not determine an optimal combination.")

if __name__ == '__main__':
    find_optimal_catalyst()