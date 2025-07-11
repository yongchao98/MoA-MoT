def recommend_htn_medications():
    """
    This function analyzes the patient's case and recommends a 3-drug regimen
    for resistant hypertension based on clinical guidelines and the patient's specific profile.
    """

    # Patient Profile Summary:
    # - Diagnosis: Resistant Stage II Hypertension (BP: 145-153/85-89)
    # - Demographics: 27-year-old African American female
    # - Comorbidities: Diabetes (A1C 6.7%), Hypercholesterolemia, Hypothyroidism, Tachycardia (HR 91)
    # - Excluded Medications: eplerenone, finerenone, HCTZ, indapamide, loop diuretics, metolazone, verapamil

    # Medication Selection Rationale:
    # 1. Lisinopril (ACE Inhibitor): Addresses the renin-angiotensin system.
    #    It is a first-line agent and is recommended for patients with diabetes for kidney protection.
    # 2. Amlodipine (Calcium Channel Blocker): A first-line agent for African American patients.
    #    It is a potent vasodilator.
    # 3. Spironolactone (Mineralocorticoid Receptor Antagonist): The recommended add-on therapy for
    #    resistant hypertension. It acts as a potassium-sparing diuretic and targets aldosterone,
    #    a key factor in treatment-resistant cases.

    medication_1 = "Lisinopril"
    medication_2 = "Amlodipine"
    medication_3 = "Spironolactone"

    print("Based on the patient's profile of resistant hypertension, the following 3-drug regimen is recommended to maximize treatment:")
    print(f"1. {medication_1}")
    print(f"2. {medication_2}")
    print(f"3. {medication_3}")

recommend_htn_medications()
<<<Lisinopril, Amlodipine, Spironolactone>>>