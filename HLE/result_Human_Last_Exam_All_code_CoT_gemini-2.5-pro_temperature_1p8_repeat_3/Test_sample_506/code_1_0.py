def recommend_htn_medications():
    """
    This script recommends and prints a 3-drug regimen for the patient's resistant hypertension
    based on their clinical profile and medication restrictions.
    """
    patient_info = {
        "current_bp": "145-153/85-89 mmHg",
        "goal_bp": "<130/90 mmHg",
        "a1c": "6.7%",
        "relevant_contraindications": "multiple diuretics (thiazide, loop) and some MRAs."
    }

    recommendations = [
        {
            "drug_name": "Valsartan",
            "drug_class": "Angiotensin II Receptor Blocker (ARB)",
            "rationale": (f"An ARB is recommended due to the patient's compelling indication of diabetes (A1C = {patient_info['a1c']}) "
                          "for renal protection and blood pressure management.")
        },
        {
            "drug_name": "Amlodipine",
            "drug_class": "Dihydropyridine Calcium Channel Blocker (CCB)",
            "rationale": ("A CCB is a first-line agent highly effective in African American patients. "
                          "This provides powerful vasodilation to help lower the current BP of {0} "
                          "towards the goal of {1}.".format(patient_info['current_bp'], patient_info['goal_bp']))
        },
        {
            "drug_name": "Spironolactone",
            "drug_class": "Mineralocorticoid Receptor Antagonist (MRA) / Potassium-Sparing Diuretic",
            "rationale": ("The standard of care for resistant hypertension, especially when other primary diuretics are contraindicated. "
                          "This adds a necessary third mechanism of action to maximize treatment effectiveness.")
        }
    ]

    print("Recommended 3-Medication Regimen for Resistant Hypertension:")
    print("-" * 60)
    for rec in recommendations:
        print(f"Medication: {rec['drug_name']} ({rec['drug_class']})")
        print(f"Rationale: {rec['rationale']}\n")

if __name__ == "__main__":
    recommend_htn_medications()