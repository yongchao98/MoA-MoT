import sys

def solve_neuroscience_question():
    """
    This script analyzes the functional profiles of dmPFC sub-regions to identify
    the one most selectively involved in emotion processing.
    """
    # Step 1: Represent the functional profiles of dmPFC sub-regions based on
    # established findings from neuroimaging meta-analyses. The 'primary_functions'
    # list indicates the core roles of each region.
    dmpfc_parcellations = {
        'rostroventral': {
            'primary_functions': ['emotion', 'self-referential thought', 'value'],
            'notes': 'Strongly associated with affective processing and is a key node in the brain\'s default mode network.'
        },
        'rostrodorsal': {
            'primary_functions': ['cognitive control', 'social cognition', 'decision making'],
            'notes': 'Associated with higher-level cognitive functions, including reasoning about others, but less purely affective.'
        },
        'caudal-left': {
            'primary_functions': ['motor planning', 'action', 'attention'],
            'notes': 'Functionally part of the pre-Supplementary Motor Area (pre-SMA), primarily involved in action and attention.'
        },
        'caudal-right': {
            'primary_functions': ['motor planning', 'action', 'attention'],
            'notes': 'Similar to its left counterpart, this region is mainly involved in motor-related processes.'
        },
        'occipital': {
            'primary_functions': ['vision'],
            'notes': 'A distinct brain lobe responsible for visual processing, not a parcellation of the dmPFC.'
        }
    }

    # Step 2 & 3: Analyze the data to find the region most purely activated by emotion.
    # "Purely" implies that emotion is a central function, not secondary to motor or
    # executive control functions.
    print("Analyzing functional profiles of dmPFC sub-regions...\n")
    
    target_function = 'emotion'
    best_candidate = None
    
    # We will score each region based on our criteria.
    # A region gets a high score if 'emotion' is a primary function and it is not
    # dominated by motor or cognitive control.
    scores = {}
    
    for region, data in dmpfc_parcellations.items():
        score = 0
        functions = data.get('primary_functions', [])
        
        # Criterion 1: Is 'emotion' a primary function?
        if target_function in functions:
            score += 10
            
        # Criterion 2: Is the region NOT dominated by motor functions?
        if 'motor planning' not in functions and 'action' not in functions:
            score += 5
            
        # Criterion 3: Is the region NOT dominated by pure cognitive control?
        if 'cognitive control' not in functions:
            score += 2
            
        scores[region] = score

    # Find the region with the highest score, excluding the distractor 'occipital'
    relevant_regions = {r: s for r, s in scores.items() if r != 'occipital'}
    best_candidate = max(relevant_regions, key=relevant_regions.get)

    # Step 4: Print the analysis and conclusion
    print("--- Functional Profile Analysis ---")
    for region in ['rostroventral', 'rostrodorsal', 'caudal-left', 'caudal-right']:
        print(f"Region: {region}")
        print(f"  - Primary Functions: {', '.join(dmpfc_parcellations[region]['primary_functions'])}")
        print(f"  - Assessment Score for 'Pure Emotion' processing: {scores[region]}")

    print("\n--- Conclusion ---")
    print(f"The analysis identifies the '{best_candidate}' sub-region as the best fit.")
    print("While other regions might be involved, the rostroventral dmPFC's core functional profile is centered on emotion and self-referential value, making it the most 'purely' activated by emotion processing compared to the cognitive and motor roles of the other sub-regions.")

solve_neuroscience_question()
<<<B>>>