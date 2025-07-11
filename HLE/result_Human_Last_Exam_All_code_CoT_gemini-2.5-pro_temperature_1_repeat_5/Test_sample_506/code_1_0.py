def recommend_htn_medications():
    """
    Analyzes the patient case and recommends a 3-drug regimen for resistant hypertension.
    """
    # Patient Profile Summary
    patient_ethnicity = "African American"
    patient_conditions = ["Stage II Resistant Hypertension", "Diabetes (A1C 6.7%)", "Hypercholesterolemia"]
    bp_goal = "<130/90 mmHg"
    unavailable_meds = [
        "eplerenone", "finerenone", "hydrochlorothiazide", "indapamide", 
        "bumetanide", "furosemide", "torsemide", "metolazone", "verapamil"
    ]

    # --- Rationale Development ---
    # Guideline-based approach for resistant HTN in an African American patient with diabetes,
    # considering the extensive list of unavailable medications.

    # Medication 1: Calcium Channel Blocker (CCB)
    med_1 = "Amlodipine"
    rationale_1 = ("A dihydropyridine Calcium Channel Blocker (CCB) is a first-line agent for hypertension in African "
                   "American patients. Amlodipine is a potent, long-acting, and well-tolerated choice.")

    # Medication 2: Angiotensin II Receptor Blocker (ARB)
    med_2 = "Losartan"
    rationale_2 = ("An Angiotensin II Receptor Blocker (ARB) is a preferred agent for patients with diabetes "
                   "due to its renal-protective benefits. It works synergistically with a CCB to lower blood pressure.")

    # Medication 3: Mineralocorticoid Receptor Antagonist (MRA)
    med_3 = "Spironolactone"
    rationale_3 = ("The patient has resistant hypertension and cannot take standard thiazide or loop diuretics. "
                   "Spironolactone is the guideline-recommended agent for resistant hypertension. It acts as a mild "
                   "diuretic and is highly effective in this scenario, and it is not on the patient's exclusion list.")

    # --- Final Output ---
    print("### Recommended 3-Drug Regimen for Resistant Hypertension ###\n")
    print("Based on the patient's profile and medication constraints, the following regimen is recommended to maximize blood pressure control:\n")

    print(f"1. Medication: {med_1}")
    print(f"   Rationale: {rationale_1}\n")

    print(f"2. Medication: {med_2}")
    print(f"   Rationale: {rationale_2}\n")

    print(f"3. Medication: {med_3}")
    print(f"   Rationale: {rationale_3}\n")

# Execute the function to print the recommendation
recommend_htn_medications()