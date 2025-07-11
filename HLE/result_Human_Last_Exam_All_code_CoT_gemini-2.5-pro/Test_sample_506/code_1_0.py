def recommend_htn_medications():
    """
    This function recommends a three-drug regimen for the patient MM's
    resistant hypertension based on her clinical profile and restrictions.
    """

    # Medication choices based on clinical guidelines and patient specifics
    medication_1 = {
        "name": "Amlodipine",
        "class": "Dihydropyridine Calcium Channel Blocker",
        "rationale": "A first-line agent for hypertension in African American patients, acting as a potent vasodilator."
    }

    medication_2 = {
        "name": "Lisinopril",
        "class": "ACE Inhibitor",
        "rationale": "Provides essential renal protection for patients with diabetes and is a cornerstone of hypertension therapy."
    }

    medication_3 = {
        "name": "Spironolactone",
        "class": "Mineralocorticoid Receptor Antagonist / Potassium-Sparing Diuretic",
        "rationale": "A guideline-recommended medication for resistant hypertension, especially when other diuretics are contraindicated."
    }
    
    recommendations = [medication_1, medication_2, medication_3]

    print("Recommended 3-Drug Regimen for Resistant Hypertension:")
    print("-" * 55)
    
    for i, med in enumerate(recommendations, 1):
        print(f"Recommendation {i}:")
        print(f"  Medication: {med['name']}")
        print(f"  Rationale: {med['rationale']}")
        print("-" * 55)

if __name__ == "__main__":
    recommend_htn_medications()