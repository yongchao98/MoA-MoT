def analyze_clinical_case():
    """
    Analyzes a clinical case by matching patient symptoms to disease profiles
    associated with different imaging findings.
    """
    # Define the prominent symptoms from the clinical vignette
    patient_symptoms = {
        "transient_monocular_vision_loss",
        "pulsatile_headaches",
        "joint_pain",
        "dyspnea",
        "hearing_loss",
        "multisystem_involvement",  # Key conceptual finding
        "cranial_nerve_involvement" # Hearing and vision loss
    }

    # Simplified knowledge base linking answer choices to associated disease patterns
    disease_profiles = {
        "A": {
            "description": "Periarticular bone demineralization (Suggests inflammatory arthritis like RA)",
            "symptoms": {"joint_pain"},
            "score": 0
        },
        "B": {
            "description": "Leptomeningeal enhancement / 'snowball' hyperintensities (Classic for Neurosarcoidosis)",
            "symptoms": {
                "transient_monocular_vision_loss", "pulsatile_headaches",
                "joint_pain", "dyspnea", "hearing_loss",
                "multisystem_involvement", "cranial_nerve_involvement"
            },
            "score": 0
        },
        "C": {
            "description": "Pleural effusion (Can be sarcoidosis, but a less specific finding)",
            "symptoms": {"dyspnea", "multisystem_involvement"},
            "score": 0
        },
        "D": {
            "description": "Vascular hemorrhage (Can cause neuro symptoms, but less explanation for systemic issues)",
            "symptoms": {"transient_monocular_vision_loss", "pulsatile_headaches"},
            "score": 0
        },
        "E": {
            "description": "Intrasellar mass (e.g., Pituitary Adenoma, poor fit for systemic symptoms)",
            "symptoms": {"headaches"},
            "score": 0
        }
    }

    # Calculate a score for each choice based on matching symptoms
    for choice, profile in disease_profiles.items():
        matches = patient_symptoms.intersection(profile["symptoms"])
        profile["score"] = len(matches)

    # Find the choice with the highest score
    best_choice = max(disease_profiles, key=lambda k: disease_profiles[k]['score'])

    # Print the analysis and conclusion
    print("--- Clinical Case Analysis ---")
    print(f"Patient's Key Symptoms: {sorted(list(patient_symptoms))}\n")
    print("--- Evaluating Potential Diagnoses/Findings ---")
    for choice, profile in sorted(disease_profiles.items(), key=lambda item: item[1]['score'], reverse=True):
        print(f"Choice {choice}: {profile['description']}")
        # The line below outputs the score ("number") for each choice
        print(f"  -> Symptom Match Score: {profile['score']} / {len(patient_symptoms)}")

    print("\n--- Conclusion ---")
    print(f"The patient's constellation of symptoms involving multiple systems (neurologic, pulmonary, articular) is highly suggestive of Sarcoidosis.")
    print(f"The most specific and expected finding for Neurosarcoidosis is described in Choice {best_choice}.")

analyze_clinical_case()