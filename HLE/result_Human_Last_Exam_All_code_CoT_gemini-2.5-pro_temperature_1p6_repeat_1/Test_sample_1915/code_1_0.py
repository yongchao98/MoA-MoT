def solve_medical_case():
    """
    Analyzes a clinical case to identify the damaged anatomical structure.
    """
    patient_symptoms = {
        "pupillary_light_reflex": False,
        "adduction": False,  # Movement towards the nose
        "elevation": False,  # Upward movement
        "depression": False, # Downward movement
    }

    anatomical_structures = {
        "A": {
            "name": "Cranial nerve VI",
            "functions": ["abduction"],
            "location": "Pons"
        },
        "B": {
            "name": "Cranial nerve VII",
            "functions": ["facial_expression", "lacrimation"],
            "location": "Pons"
        },
        "C": {
            "name": "Reticular formation",
            "functions": ["arousal", "consciousness"],
            "location": "Brainstem"
        },
        "D": {
            "name": "Medulla oblongata",
            "functions": ["swallowing", "respiration"],
            "contains_nuclei": ["CN IX", "CN X", "CN XI", "CN XII"]
        },
        "E": {
            "name": "Midbrain",
            "functions": ["motor_control", "auditory_processing", "visual_processing"],
            "contains_nuclei": ["CN III", "CN IV"]
        }
    }

    # Functions controlled by CN III, whose nucleus is in the Midbrain
    cn_iii_functions = {
        "adduction": True,
        "elevation": True,
        "depression": True, # Partial, via inferior rectus
        "pupillary_light_reflex": True # Efferent limb
    }

    print("Analyzing patient symptoms:")
    print(f"- Pupillary Light Reflex: {'Intact' if patient_symptoms['pupillary_light_reflex'] else 'Absent'}")
    print(f"- Adduction: {'Possible' if patient_symptoms['adduction'] else 'Not Possible'}")
    print(f"- Elevation: {'Possible' if patient_symptoms['elevation'] else 'Not Possible'}")
    print(f"- Depression: {'Possible' if patient_symptoms['depression'] else 'Not Possible'}")
    print("\nConclusion from symptoms: The pattern of deficits strongly suggests a Cranial Nerve III palsy.\n")

    print("Evaluating anatomical structures:")
    best_match = None
    max_score = 0

    for key, structure in anatomical_structures.items():
        score = 0
        reasoning = []
        if "contains_nuclei" in structure and "CN III" in structure["contains_nuclei"]:
            reasoning.append(f"Contains the nucleus for Cranial Nerve III.")
            # Check how many of the patient's symptoms are explained
            for symptom, affected in patient_symptoms.items():
                if not affected and cn_iii_functions.get(symptom):
                    score += 1
            if score > max_score:
                max_score = score
                best_match = key
        else:
            reasoning.append("Does not contain the CN III nucleus, which is responsible for the observed symptoms.")

        print(f"Choice {key}: {structure['name']}")
        print(f"  - Reasoning: {' '.join(reasoning)}")


    print("\n---\nFinal Determination:")
    if best_match:
        print(f"The structure that best explains the patient's complete CN III palsy is the '{anatomical_structures[best_match]['name']}', as it houses the CN III nucleus.")
        print(f"The correct answer is therefore {best_match}.")
    else:
        print("Could not determine the best match based on the provided information.")

solve_medical_case()