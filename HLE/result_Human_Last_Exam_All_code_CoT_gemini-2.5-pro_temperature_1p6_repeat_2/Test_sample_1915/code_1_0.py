def solve_clinical_case():
    """
    Analyzes clinical symptoms to identify the damaged anatomical structure.
    """

    # Define the patient's symptoms based on the clinical description
    patient_symptoms = {
        "Pupillary Light Reflex": "Absent (implicates efferent CN III)",
        "Adduction (inward eye movement)": "Unable (implicates CN III)",
        "Depression (downward eye movement)": "Unable (implicates CN III)",
        "Elevation (upward eye movement)": "Unable (implicates CN III)"
    }

    # Map cranial nerves to their primary functions and anatomical location of their nucleus
    anatomical_structures = {
        "A. Cranial nerve VI": {"functions": ["Abduction"], "location": "Pons"},
        "B. Cranial nerve VII": {"functions": ["Facial expression", "Taste"], "location": "Pons"},
        "C. Reticular formation": {"functions": ["Arousal", "Consciousness"], "location": "Brainstem (diffuse)"},
        "D. Medulla oblongata": {"functions": ["Swallowing", "Heart rate"], "location_of_nuclei": ["CN IX", "CN X", "CN XI", "CN XII"]},
        "E. Midbrain": {"functions": ["Eye movement", "Pupillary reflex"], "location_of_nuclei": ["CN III", "CN IV"]}
    }

    # The symptoms consistently point to a deficit in Cranial Nerve III
    implicated_nerve = "CN III"

    print("--- Clinical Reasoning Steps ---")
    print("1. Analyzing Patient's Symptoms:")
    for symptom, status in patient_symptoms.items():
        print(f"   - {symptom}: {status}")

    print(f"\n2. Synthesizing Findings:")
    print(f"   All presented symptoms are deficits related to the functions of Cranial Nerve {implicated_nerve}.")

    print("\n3. Locating the Anatomical Structure:")
    correct_answer_key = None
    correct_structure_name = None
    for key, value in anatomical_structures.items():
        if "location_of_nuclei" in value and implicated_nerve in value["location_of_nuclei"]:
            correct_answer_key = key.split('.')[0]
            correct_structure_name = key.split('. ')[1]
            break

    print(f"   The nucleus of {implicated_nerve} is located in the {correct_structure_name}.")
    print(f"   Therefore, damage to the {correct_structure_name} explains the patient's presentation.")

    print("\n--- Final Answer ---")
    print(f"The patient's presentation is explained by damage to the {correct_structure_name}.")
    # The format required is <<<AnswerKey>>>
    print(f"\nFinal Answer Choice: <<<E>>>")


solve_clinical_case()