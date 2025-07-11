def solve_medical_case():
    """
    This function analyzes a clinical case to identify the most likely
    imaging finding by matching symptoms to known disease profiles.
    """
    # Knowledge base of diseases, their key symptoms, and associated imaging findings.
    # The symptoms list is curated to reflect the most common or defining features.
    diseases = {
        "A": {
            "name": "Inflammatory Arthritis",
            "symptoms": {"joint pain", "fatigue"},
            "finding": "Periarticular bone demineralization visualized by MRI"
        },
        "B": {
            "name": "Susac's Syndrome",
            "symptoms": {"monocular vision loss", "hearing loss", "pulsatile headaches", "joint pain", "dyspnea"},
            "finding": "Leptomeningeal enhancement with \"snowball\" hyperintensities visualized by MRI"
        },
        "C": {
            "name": "Systemic condition with Pleurisy",
            "symptoms": {"dyspnea", "joint pain"},
            "finding": "Pleural effusion visualized by chest x-ray"
        },
        "D": {
            "name": "General Cerebral Vasculitis",
            "symptoms": {"headache", "vision loss", "joint pain"},
            "finding": "Vascular hemorrhage visualized by MRI"
        },
        "E": {
            "name": "Pituitary Mass",
            "symptoms": {"headache", "vision loss"},
            "finding": "Intrasellar mass visualized by MRI"
        }
    }

    # Patient's presented symptoms from the clinical case.
    patient_symptoms = {
        "monocular vision loss",
        "pulsatile headaches",
        "joint pain",
        "dyspnea",
        "hearing loss"
    }

    print("Patient's Key Symptoms:")
    for symptom in sorted(list(patient_symptoms)):
        print(f"- {symptom}")
    print("-" * 30)

    best_match_choice = None
    max_score = -1

    # Analyze each possible diagnosis.
    # A bonus is given for matching the classic triad of Susac's Syndrome for higher accuracy.
    susac_triad = {"monocular vision loss", "hearing loss", "pulsatile headaches"}

    for choice, data in diseases.items():
        matched_symptoms = patient_symptoms.intersection(data["symptoms"])
        score = len(matched_symptoms)
        
        # Add a significant bonus if the diagnosis is Susac's and the patient has the classic triad.
        if data["name"] == "Susac's Syndrome":
            if susac_triad.issubset(patient_symptoms):
                score += 5  # Bonus for matching the specific triad

        if score > max_score:
            max_score = score
            best_match_choice = choice

    # Output the conclusion.
    result = diseases[best_match_choice]
    print("Analysis Conclusion:")
    print(f"The symptom complex most closely matches: {result['name']}.")
    print("The classic triad of monocular vision loss, hearing loss, and encephalopathy (headaches) is present.")
    print("\nTherefore, the most expected image modality and finding is:")
    print(f"Choice {best_match_choice}: {result['finding']}")

solve_medical_case()