def recommend_htn_medications():
    """
    This function recommends a 3-drug regimen for the patient's resistant hypertension
    based on clinical guidelines and the provided patient information.
    """
    # Rationale:
    # 1. Chlorthalidone: A potent thiazide-like diuretic, not on the patient's exclusion list,
    #    and is a cornerstone for treating resistant hypertension.
    # 2. Amlodipine: A dihydropyridine calcium channel blocker (CCB), a first-line agent for
    #    African American patients, and works well in combination.
    # 3. Lisinopril: An ACE inhibitor, recommended for patients with hypertension and diabetes
    #    to provide renal protection.
    
    medication_1 = "Chlorthalidone"
    medication_2 = "Amlodipine"
    medication_3 = "Lisinopril"

    print("Recommended 3-drug regimen for resistant hypertension:")
    print(f"1. {medication_1}")
    print(f"2. {medication_2}")
    print(f"3. {medication_3}")

recommend_htn_medications()