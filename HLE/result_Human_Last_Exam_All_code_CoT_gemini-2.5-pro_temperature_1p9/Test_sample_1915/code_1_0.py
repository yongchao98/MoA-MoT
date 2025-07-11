def solve_medical_case():
    """
    Analyzes patient symptoms to identify the damaged anatomical structure.
    """
    # Patient's key symptoms in the right eye
    symptoms = {
        "no_pupillary_light_reflex",
        "inability_to_adduct",
        "inability_to_depress",
        "inability_to_elevate"
    }

    # Map symptoms to the controlling cranial nerve
    symptom_map = {
        "no_pupillary_light_reflex": "CN III (Parasympathetic fibers)",
        "inability_to_adduct": "CN III (Oculomotor Nerve)",
        "inability_to_depress": "CN III (Oculomotor Nerve)",
        "inability_to_elevate": "CN III (Oculomotor Nerve)"
    }
    
    # Relevant anatomical structures from the choices and their functions
    anatomical_structures = {
        "A": {"name": "Cranial nerve VI", "function": "Controls eye abduction (outward movement). Does not match symptoms."},
        "B": {"name": "Cranial nerve VII", "function": "Controls muscles of facial expression. Does not match symptoms."},
        "C": {"name": "Reticular formation", "function": "Diffuse network for arousal. Less specific than other choices."},
        "D": {"name": "Medulla oblongata", "function": "Origin of lower cranial nerves (IX, X, XI, XII). Does not match symptoms."},
        "E": {"name": "Midbrain", "function": "Contains the nucleus of Cranial Nerve III. A lesion here explains all symptoms."}
    }

    print("Step 1: Analyzing the patient's symptoms.")
    print("The patient's deficits are:")
    for symptom in symptoms:
        print(f"- {symptom.replace('_', ' ').capitalize()}")
    print("\nAll these functions are controlled by a single cranial nerve.")
    
    print("\nStep 2: Identifying the affected cranial nerve.")
    cranial_nerve_number = 3
    print(f"The symptoms collectively point to a complete palsy of Cranial Nerve {cranial_nerve_number} (the Oculomotor Nerve).")
    
    print("\nStep 3: Correlating the affected nerve to its anatomical origin from the choices.")
    print(f"We need to find the structure that contains the nucleus for Cranial Nerve {cranial_nerve_number}.")
    for choice, data in anatomical_structures.items():
        print(f"- {choice}. {data['name']}: {data['function']}")
    
    print("\nStep 4: Final Conclusion.")
    correct_choice = "E"
    conclusion = f"The Midbrain is the location of the Oculomotor Nerve (CN {cranial_nerve_number}) nucleus. Therefore, damage to the {anatomical_structures[correct_choice]['name']} best explains the patient's presentation."
    print(conclusion)

solve_medical_case()