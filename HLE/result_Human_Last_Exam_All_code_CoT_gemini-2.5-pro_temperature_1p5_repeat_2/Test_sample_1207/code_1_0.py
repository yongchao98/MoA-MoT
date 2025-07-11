def diagnose_patient_case():
    """
    Analyzes a clinical case to determine the most likely imaging finding.
    """
    patient_profile = {
        "age": 44,
        "sex": "female",
        "key_symptoms": [
            "transient monocular vision loss",
            "pulsatile headaches",
            "joint pain",
            "dyspnea",
            "hearing loss",
            "painful lower extremity lesion"
        ]
    }

    print("Step 1: Analyzing Patient Symptoms")
    print(f"Patient is a {patient_profile['age']}-year-old {patient_profile['sex']} with multi-system complaints:")
    for symptom in patient_profile['key_symptoms']:
        print(f"- {symptom}")
    print("\nThis clinical picture suggests a systemic inflammatory or granulomatous disease.\n")

    print("Step 2: Evaluating Differential Diagnoses and Imaging Findings")
    diagnoses = {
        "A": {
            "finding": "Periarticular bone demineralization by MRI",
            "disease": "Rheumatoid Arthritis",
            "fit": "Poor fit. Does not explain neurological or pulmonary symptoms well."
        },
        "B": {
            "finding": "Leptomeningeal enhancement with 'snowball' hyperintensities by MRI",
            "disease": "Neurosarcoidosis",
            "fit": "Excellent fit. This is a classic finding for the neurological manifestation of Sarcoidosis, which explains the full constellation of multi-system symptoms (ocular, neuro, pulmonary, joint, skin)."
        },
        "C": {
            "finding": "Pleural effusion by chest x-ray",
            "disease": "Multiple causes, including Sarcoidosis",
            "fit": "Possible, but not the most specific finding and doesn't address the chief neurological complaints."
        },
        "D": {
            "finding": "Vascular hemorrhage by MRI",
            "disease": "Vasculitis, hypertension",
            "fit": "Possible, but inflammatory changes are more characteristic for the suspected underlying disease process than hemorrhage."
        },
        "E": {
            "finding": "Intrasellar mass by MRI",
            "disease": "Pituitary Adenoma",
            "fit": "Poor fit. Cannot explain the systemic symptoms like dyspnea, joint pain, and hearing loss."
        }
    }

    best_fit_choice = ""
    best_fit_reason = ""
    max_fit_score = 0

    print("Comparing answer choices to the clinical picture:")
    for choice, data in diagnoses.items():
        print(f"Choice {choice}: {data['finding']}")
        print(f"  Associated with: {data['disease']}")
        print(f"  Assessment: {data['fit']}\n")
        # A simple way to score the fit for demonstration
        if "Excellent fit" in data['fit']:
            score = 3
        elif "Possible" in data['fit']:
            score = 2
        else:
            score = 1
        
        if score > max_fit_score:
            max_fit_score = score
            best_fit_choice = choice
            best_fit_reason = data['fit']

    print("Step 3: Final Conclusion")
    print("The patient's presentation of multi-system inflammation strongly points to Sarcoidosis.")
    print("The most specific and expected imaging finding for the neurological symptoms described is characteristic of Neurosarcoidosis.")
    print(f"\nThe best matching option is {best_fit_choice} because: {best_fit_reason}")

diagnose_patient_case()