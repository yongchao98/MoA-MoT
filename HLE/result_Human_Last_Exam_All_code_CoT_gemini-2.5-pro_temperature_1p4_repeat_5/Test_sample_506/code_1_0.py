def recommend_htn_medication_regimen():
    """
    Analyzes the patient case and recommends a 3-drug regimen for resistant hypertension.
    The patient's intolerances are considered in this recommendation.
    """
    # Patient's key clinical numbers
    patient_bp_systolic_max = 153
    patient_bp_diastolic_max = 89
    goal_bp_systolic = 130
    goal_bp_diastolic = 90

    # Medication Selection based on clinical guidelines and patient profile
    
    # 1. Dihydropyridine Calcium Channel Blocker (CCB)
    # Rationale: First-line agent for African American patients. Not on the exclusion list.
    medication_1 = "Amlodipine"
    
    # 2. Angiotensin II Receptor Blocker (ARB)
    # Rationale: Preferred for patients with diabetes for renal protection. Well-tolerated.
    medication_2 = "Valsartan"
    
    # 3. Mineralocorticoid Receptor Antagonist (MRA)
    # Rationale: Guideline-recommended agent for resistant hypertension. Spironolactone is not on the exclusion list.
    medication_3 = "Spironolactone"

    print("Recommendation for Resistant Hypertension Treatment Maximization:")
    print("---------------------------------------------------------------")
    print(f"The patient's current blood pressure is as high as {patient_bp_systolic_max}/{patient_bp_diastolic_max} mmHg, with a goal of <{goal_bp_systolic}/{goal_bp_diastolic} mmHg.")
    print("Given the patient's profile and medication constraints, the following three medications are recommended to form an optimal regimen:\n")
    
    print(f"1. {medication_1} (Calcium Channel Blocker)")
    print(f"2. {medication_2} (Angiotensin II Receptor Blocker)")
    print(f"3. {medication_3} (Mineralocorticoid Receptor Antagonist/Diuretic)")

    print("\nThis combination targets three different pathophysiological pathways and is a standard-of-care approach for resistant hypertension.")
    
    # This line formats the patient's BP numbers and recommended drugs into a conceptual equation as requested.
    print(f"\nConceptual Goal Equation: ({patient_bp_systolic_max}/{patient_bp_diastolic_max}) + [{medication_1} + {medication_2} + {medication_3}] -> Goal BP (<{goal_bp_systolic}/{goal_bp_diastolic})")


# Execute the function to display the recommendations
recommend_htn_medication_regimen()