def recommend_htn_medications():
    """
    This function outlines and prints a recommended 3-drug regimen for a patient 
    with resistant hypertension, considering their specific clinical profile and constraints.
    """
    
    patient_name = "MM"
    condition = "Resistant Stage II Hypertension with Type 2 Diabetes"

    print(f"Recommended 3-drug hypertension regimen for patient {patient_name} with {condition}:\n")

    # Medication 1: Calcium Channel Blocker
    med_1_name = "Amlodipine"
    med_1_class = "Dihydropyridine Calcium Channel Blocker (CCB)"
    med_1_rationale = ("This is a first-line agent for African American patients. "
                     "It is a potent vasodilator effective for lowering blood pressure.")

    # Medication 2: Angiotensin II Receptor Blocker
    med_2_name = "Losartan"
    med_2_class = "Angiotensin II Receptor Blocker (ARB)"
    med_2_rationale = ("Indicated due to the patient's diabetes for renal protection. "
                       "It blocks the RAAS pathway, a key BP regulator.")

    # Medication 3: Mineralocorticoid Receptor Antagonist
    med_3_name = "Spironolactone"
    med_3_class = "Mineralocorticoid Receptor Antagonist (MRA) / Potassium-Sparing Diuretic"
    med_3_rationale = ("Guideline-recommended for resistant hypertension. It provides a diuretic "
                       "effect, which is crucial since other standard diuretics are not an option for this patient.")

    print("1. Medication: " + med_1_name)
    print("   Class: " + med_1_class)
    print("   Rationale: " + med_1_rationale)
    print("-" * 70)

    print("2. Medication: " + med_2_name)
    print("   Class: " + med_2_class)
    print("   Rationale: " + med_2_rationale)
    print("-" * 70)
    
    print("3. Medication: " + med_3_name)
    print("   Class: " + med_3_class)
    print("   Rationale: " + med_3_rationale)
    print("-" * 70)

recommend_htn_medications()