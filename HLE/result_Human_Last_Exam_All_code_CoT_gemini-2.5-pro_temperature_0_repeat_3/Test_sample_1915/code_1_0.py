def diagnose_neurological_lesion():
    """
    Analyzes patient symptoms to identify the damaged anatomical structure.
    """
    # Step 1: Define the patient's key symptoms for the right eye.
    patient_symptoms = {
        "no pupillary light reflex",
        "inability to adduct",
        "inability to depress",
        "inability to elevate"
    }

    # Step 2 & 3: Define the functions/symptoms associated with each anatomical structure.
    # This simplified model maps structures to the symptoms their damage would cause.
    anatomical_structures = {
        "A. Cranial nerve VI": {
            "description": "Innervates the lateral rectus muscle, responsible for eye abduction (moving eye outwards).",
            "symptoms": {"inability to abduct"}
        },
        "B. Cranial nerve VII": {
            "description": "The Facial Nerve controls muscles of facial expression and taste.",
            "symptoms": {"facial paralysis", "loss of taste"}
        },
        "C. Reticular formation": {
            "description": "A complex network in the brainstem involved in consciousness and arousal.",
            "symptoms": {"altered consciousness", "coma"}
        },
        "D. Medulla oblongata": {
            "description": "The lowest part of the brainstem, containing nuclei for CN IX, X, XI, XII.",
            "symptoms": {"difficulty swallowing", "hoarseness", "irregular breathing"}
        },
        "E. Midbrain": {
            "description": "The upper part of the brainstem containing the nucleus for Cranial Nerve III (Oculomotor).",
            "symptoms": {
                "no pupillary light reflex",  # CN III parasympathetic function
                "inability to elevate",      # CN III motor function (superior rectus, inferior oblique)
                "inability to depress",      # CN III motor function (inferior rectus)
                "inability to adduct"       # CN III motor function (medial rectus)
            }
        }
    }

    # Step 4: Analyze which structure best explains the symptoms.
    best_match = None
    max_match_count = -1

    print("--- Patient Symptom Analysis ---")
    print(f"Patient presents with the following symptoms in the right eye: {', '.join(sorted(list(patient_symptoms)))}.")
    print("These symptoms point to a complete palsy of the Oculomotor Nerve (Cranial Nerve III).\n")
    print("--- Evaluating Potential Lesion Locations ---")

    for structure_key, data in anatomical_structures.items():
        # Check how many of the patient's symptoms are explained by damage to this structure.
        explained_symptoms = patient_symptoms.intersection(data["symptoms"])
        match_count = len(explained_symptoms)

        print(f"Location: {structure_key}")
        print(f"  - Function/Relevance: {data['description']}")
        print(f"  - Explains {match_count} of {len(patient_symptoms)} patient symptoms.")
        if match_count > max_match_count:
            max_match_count = match_count
            best_match = structure_key

    # Step 5: Print the conclusion.
    print("\n--- Conclusion ---")
    print(f"The location that fully explains all the patient's symptoms is the '{best_match}'.")
    print("Damage to the Midbrain is the only option that accounts for a complete Cranial Nerve III palsy, which includes loss of pupillary reflex and the inability to adduct, elevate, and depress the eye.")

diagnose_neurological_lesion()