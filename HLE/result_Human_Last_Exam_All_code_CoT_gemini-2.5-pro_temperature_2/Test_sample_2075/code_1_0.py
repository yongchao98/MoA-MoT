import json

def find_most_specialized_region():
    """
    This function identifies the dmPFC sub-region most purely specialized
    for emotion processing based on established neuroscientific findings.
    """
    # This data structure summarizes findings from functional neuroimaging
    # meta-analyses on the specialization of medial prefrontal cortex sub-regions.
    # The 'primary_function' field describes the core role of each area.
    dmpfc_functions = {
        "caudal-right": {
            "primary_function": "Cognitive Control & Action Selection",
            "notes": "Closely tied to motor planning and execution. Activated by tasks requiring response selection and inhibition."
        },
        "rostroventral": {
            "primary_function": "Emotion & Valuation",
            "notes": "Core region for processing affect, assigning emotional value to stimuli, and autonomic regulation. Tightly linked to limbic structures."
        },
        "rostrodorsal": {
            "primary_function": "Social Cognition & Mentalizing",
            "notes": "Involved in thinking about the mental states of others (Theory of Mind) and self-referential thought."
        },
        "caudal-left": {
            "primary_function": "Cognitive Control & Action Selection",
            "notes": "Similar to its right hemisphere counterpart, involved in high-level control of action."
        },
        "occipital": {
            "primary_function": "Vision",
            "notes": "This is a different brain lobe and not part of the dmPFC. It is the primary visual processing center."
        }
    }

    target_function = "Emotion"
    most_purely_emotional_region = None
    justification = ""

    # We search for the region whose primary function is explicitly 'Emotion & Valuation'.
    # This aligns with the 'most purely activated' criterion, as it distinguishes
    # core emotional processing from the cognitive or action-related components of
    # emotional experiences that might activate other regions.
    for region, data in dmpfc_functions.items():
        if target_function in data["primary_function"]:
            most_purely_emotional_region = region
            justification = data["notes"]
            break

    if most_purely_emotional_region:
        print(f"Based on functional specialization, the sub-region most purely activated by emotion processing is:")
        print(f"\nRegion: {most_purely_emotional_region}")
        print(f"Primary Function: {dmpfc_functions[most_purely_emotional_region]['primary_function']}")
        print(f"Justification: {justification}")
    else:
        print("Could not identify a region specialized for the target function.")

find_most_specialized_region()