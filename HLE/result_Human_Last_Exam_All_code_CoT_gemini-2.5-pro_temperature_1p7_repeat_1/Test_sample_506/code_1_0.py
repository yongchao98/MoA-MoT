def recommend_htn_medications():
    """
    Analyzes the patient case and recommends three medications for resistant hypertension,
    considering the patient's profile and medication constraints.
    """

    patient_profile = {
        "Ethnicity": "African American",
        "Comorbidities": ["Type 2 Diabetes", "Hypercholesterolemia", "Resistant Stage II HTN"],
        "BP_Goal": "<130/90 mmHg",
        "Excluded_Medications": [
            "eplerenone", "finerenone", "hydrochlorothiazide", "indapamide",
            "bumetanide", "furosemide", "torsemide", "metolazone", "verapamil"
        ]
    }

    # Based on the plan, select three optimal medications.
    recommendations = {
        "Losartan": "An Angiotensin II Receptor Blocker (ARB). This is recommended due to the patient's diabetes for kidney protection and blood pressure control.",
        "Amlodipine": "A dihydropyridine Calcium Channel Blocker (CCB). This is a first-line agent for African American patients and is a potent vasodilator.",
        "Chlorthalidone": "A long-acting thiazide-like diuretic. This drug is not on the patient's exclusion list and is highly effective for resistant hypertension by reducing blood volume."
    }

    print("Recommended 3-Drug Regimen for Maximal Hypertension Treatment:")
    print("-" * 60)
    for med, reason in recommendations.items():
        print(f"1. Medication: {med}")
        print(f"   Rationale: {reason}\n")

recommend_htn_medications()