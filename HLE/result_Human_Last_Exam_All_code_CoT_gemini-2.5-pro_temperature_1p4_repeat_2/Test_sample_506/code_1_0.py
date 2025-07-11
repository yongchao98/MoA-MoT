import pandas as pd

def recommend_htn_medications():
    """
    Recommends a 3-drug regimen for the patient's resistant hypertension
    based on clinical profile and medication exclusions.
    """
    # Patient Profile
    # Race: African American
    # Comorbidities: Type 2 Diabetes, Hypercholesterolemia
    # Condition: Stage II Resistant Hypertension
    # Excluded Medications: eplerenone, finerenone, HCTZ, indapamide, loop diuretics, verapamil

    # Drug Selection Rationale:
    # 1. ARB (Losartan): Recommended for patients with diabetes for renal protection.
    # 2. Dihydropyridine CCB (Amlodipine): First-line for African American patients.
    # 3. MRA/Diuretic (Spironolactone): Guideline-recommended for resistant hypertension,
    #    especially when other diuretics are not an option. It is not on the excluded list.

    recommendations = {
        "Medication": ["Losartan", "Amlodipine", "Spironolactone"],
        "Drug Class": ["Angiotensin II Receptor Blocker (ARB)", "Dihydropyridine Calcium Channel Blocker (CCB)", "Mineralocorticoid Receptor Antagonist (MRA) / Diuretic"],
        "Rationale": [
            "Addresses renin-angiotensin system; provides kidney protection for diabetes.",
            "First-line agent for African American patients; potent vasodilator.",
            "Effective for resistant hypertension; circumvents the specified diuretic exclusions."
        ]
    }

    recommendation_df = pd.DataFrame(recommendations)

    print("Recommended 3-Drug Regimen for Hypertension:\n")
    print(recommendation_df.to_string(index=False))

if __name__ == '__main__':
    recommend_htn_medications()