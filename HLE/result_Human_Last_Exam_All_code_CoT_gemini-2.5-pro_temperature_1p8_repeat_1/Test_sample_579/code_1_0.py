import pandas as pd

def analyze_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    The function uses a simple scoring model based on key features.
    """

    # Patient data extracted from the vignette
    patient_data = {
        "age": 64,
        "bmi": 39,
        "is_smoker": True,
        "smoking_years": 15,
        "lesion_locations": ["axillary", "inframammary", "inguinal"],
        "lesion_morphology": ["bullae", "plaques", "purulent nodules"]
    }

    # Criteria for each diagnosis
    diagnoses = {
        "A. Malignant Intertrigo": {
            "locations": ["inframammary", "axillary"],
            "morphology": ["plaques", "nodules"],
            "risk_factors": ["history_of_cancer"],
            "contraindications": ["purulent nodules", "bullae"]
        },
        "B. Allergic contact dermatitis": {
            "locations": ["any"],
            "morphology": ["plaques", "bullae", "vesicles"],
            "risk_factors": [],
            "contraindications": ["purulent nodules"]
        },
        "C. Hidradenitis Suppurativa": {
            "locations": ["axillary", "inframammary", "inguinal"],
            "morphology": ["purulent nodules", "plaques", "abscesses"],
            "risk_factors": ["obesity", "smoking"],
            "contraindications": []
        },
        "D. Atopic dermatitis": {
            "locations": ["flexural"],
            "morphology": ["plaques", "lichenification"],
            "risk_factors": [],
            "contraindications": ["purulent nodules", "bullae"]
        },
        "E. Psoriasis": {
            "locations": ["axillary", "inframammary", "inguinal"],
            "morphology": ["plaques"],
            "risk_factors": [],
            "contraindications": ["purulent nodules", "bullae"]
        }
    }

    scores = {}

    print("--- Diagnostic Analysis ---")

    for name, criteria in diagnoses.items():
        score = 0
        # Score based on locations
        if any(loc in patient_data["lesion_locations"] for loc in criteria["locations"]) or "any" in criteria["locations"]:
            score += 2
        if all(loc in patient_data["lesion_locations"] for loc in criteria["locations"]):
            score += 1 # Bonus for perfect match
        if "flexural" in criteria["locations"]: # Generic term for these locations
             score += 2

        # Score based on morphology
        if any(morph in patient_data["lesion_morphology"] for morph in criteria["morphology"]):
            score += 2
        if "purulent nodules" in criteria["morphology"] and "purulent nodules" in patient_data["lesion_morphology"]:
            score += 3 # Highly specific finding

        # Score based on risk factors
        if "obesity" in criteria["risk_factors"] and patient_data["bmi"] > 30:
            score += 2
        if "smoking" in criteria["risk_factors"] and patient_data["is_smoker"]:
            score += 2

        # Penalize for contraindications
        if any(morph in patient_data["lesion_morphology"] for morph in criteria["contraindications"]):
            score -= 3

        scores[name] = score

    # Determine the highest score
    best_diagnosis = max(scores, key=scores.get)

    print("\nPatient's Key Information:")
    print(f"- Age: {patient_data['age']}")
    print(f"- BMI: {patient_data['bmi']} (Class II Obesity)")
    print(f"- Smoker: Yes, for {patient_data['smoking_years']} years")
    print(f"- Lesion Locations: {', '.join(patient_data['lesion_locations'])}")
    print(f"- Lesion Types: {', '.join(patient_data['lesion_morphology'])}")

    print("\nDiagnostic Scores:")
    for name, score in scores.items():
        print(f"- {name}: {score}")

    print("\n--- Conclusion ---")
    print(f"The most likely diagnosis is {best_diagnosis}.")
    print("\nReasoning:")
    print("This diagnosis is strongly supported by the classic triad of factors present in the patient:")
    print("1. **Location:** Lesions are in the characteristic intertriginous areas (axillary, inframammary, inguinal).")
    print("2. **Morphology:** The presence of 'purulent nodules' is a hallmark sign of Hidradenitis Suppurativa.")
    print(f"3. **Risk Factors:** The patient has major risk factors, including a BMI of {patient_data['bmi']} (obesity) and a {patient_data['smoking_years']}-year smoking history.")
    print("While other conditions can affect these areas, none align as closely with the complete clinical picture, especially the purulent nodules.")

# Execute the analysis
analyze_clinical_case()