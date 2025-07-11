import collections

def analyze_clinical_case():
    """
    Analyzes a clinical presentation to identify the damaged anatomical structure.
    """
    # Patient's key symptoms related to the right eye
    symptoms = {
        "pupillary light reflex": "absent",
        "adduction": "impaired",
        "depression": "impaired",
        "elevation": "impaired"
    }

    # Map of anatomical structures to their relevant functions
    anatomical_functions = {
        "A. Cranial nerve VI": {
            "name": "Abducens Nerve",
            "controls": ["abduction"]
        },
        "B. Cranial nerve VII": {
            "name": "Facial Nerve",
            "controls": ["facial expression", "taste", "eyelid closure"]
        },
        "C. Reticular formation": {
            "name": "Reticular Formation",
            "controls": ["arousal", "consciousness", "sleep-wake cycle"]
        },
        "D. Medulla oblongata": {
            "name": "Medulla Oblongata",
            "controls": ["contains nuclei for CN IX, X, XI, XII"]
        },
        "E. Midbrain": {
            "name": "Midbrain",
            "controls": [
                "contains nucleus for CN III (Oculomotor Nerve)",
                "adduction", 
                "depression", 
                "elevation",
                "pupillary light reflex"
            ]
        }
    }
    
    # The functions lost by the patient correspond to CN III (Oculomotor nerve) functions
    lost_functions = [
        "adduction", 
        "depression", 
        "elevation", 
        "pupillary light reflex"
    ]
    
    print("Step 1: Identify the patient's neurological deficits.")
    print(f"The patient has deficits in the following functions of the right eye: {', '.join(lost_functions)}.\n")

    print("Step 2: Correlate deficits with a cranial nerve.")
    print("These functions are all controlled by the Oculomotor Nerve (Cranial Nerve III).\n")

    print("Step 3: Identify the anatomical location of the responsible nerve's nucleus.")
    
    best_match = None
    for choice, data in anatomical_functions.items():
        # Check if the structure's known functions cover all the patient's deficits
        is_match = all(func in data["controls"] for func in lost_functions)
        if is_match:
            best_match = choice
            print(f"Analyzing option {choice} ({data['name']})...")
            print(f"The {data['name']} controls the functions of the Oculomotor Nerve (CN III), which accounts for all the patient's symptoms.")
            break
            
    if best_match:
        print("\nConclusion: The patient's presentation is explained by damage to the anatomical structure that contains the Oculomotor Nerve (CN III) nucleus.")
        print(f"The correct option is {best_match.split('.')[0]}.")
    else:
        print("Could not find a direct match in the provided options.")

analyze_clinical_case()