import sys

def solve_neuro_question():
    """
    This function identifies the dmPFC sub-region most purely activated
    by emotion processing based on established neuroscience findings.
    """
    # Step 1: Define the primary functions of each brain region based on
    # large-scale meta-analyses in neuroscience. The first function in each list
    # represents the most specialized or "pure" role for that region.
    parcellation_functions = {
        'caudal-right': ['Action Selection', 'Motor Control', 'Conflict Monitoring'],
        'rostroventral': ['Emotion Processing', 'Valuation', 'Self-Referential Thought'],
        'rostrodorsal': ['Social Cognition (Theory of Mind)', 'Metacognition', 'Reasoning'],
        'caudal-left': ['Action Selection', 'Motor Control', 'Conflict Monitoring'],
        'occipital': ['Visual Processing'] # Included as a distractor
    }

    target_function = 'Emotion Processing'
    most_purely_activated_region = None
    reasoning = ""

    print("Analyzing functional specializations of dmPFC sub-regions...\n")

    # Step 2: Iterate through the regions to find the one where the target function is primary.
    for region, functions in parcellation_functions.items():
        primary_function = functions[0]
        print(f"Region: {region}")
        print(f"  - Primary Function: {primary_function}")
        print(f"  - Other Associated Functions: {', '.join(functions[1:]) if len(functions) > 1 else 'None'}\n")
        
        if primary_function.lower() == target_function.lower():
            most_purely_activated_region = region
            reasoning = (
                f"The '{region}' sub-region is identified as the area most purely "
                f"activated by '{target_function}'.\nWhile other regions might be involved in emotion, "
                f"this sub-region's primary role is dedicated to processing affective and emotional states."
            )

    # Step 3: Output the final conclusion.
    print("---Conclusion---")
    if most_purely_activated_region:
        print(reasoning)
        print(f"\nThe correct answer is: {most_purely_activated_region}")
    else:
        print(f"Could not find a region primarily associated with '{target_function}'.")

# Execute the function to find the answer.
solve_neuro_question()