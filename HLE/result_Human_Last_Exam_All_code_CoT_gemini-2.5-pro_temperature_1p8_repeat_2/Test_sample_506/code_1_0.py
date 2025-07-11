def recommend_htn_medications():
    """
    This function recommends a 3-drug regimen for the patient's resistant hypertension
    based on her clinical profile and medication restrictions.
    """
    # Rationale:
    # 1. Lisinopril (ACE Inhibitor): Standard of care for patients with hypertension and diabetes.
    # 2. Amlodipine (Calcium Channel Blocker): First-line agent, highly effective in African American patients.
    # 3. Spironolactone (Mineralocorticoid Receptor Antagonist): Guideline-recommended agent for resistant hypertension,
    #    and it is not on the patient's exclusion list.

    medication_1 = "Lisinopril"
    medication_2 = "Amlodipine"
    medication_3 = "Spironolactone"

    print("Recommended 3-Drug Regimen for Maximal Hypertension Treatment:")
    print(f"1. {medication_1}")
    print(f"2. {medication_2}")
    print(f"3. {medication_3}")

recommend_htn_medications()