def recommend_hypertension_medications():
    """
    This script analyzes the patient's case and recommends a three-drug regimen
    for resistant hypertension based on clinical guidelines and patient-specific factors.
    """

    # Patient details for context in the output
    patient_condition = "Resistant Stage II Hypertension"
    patient_comorbidities = ["Diabetes (A1C 6.7%)", "Hypercholesterolemia", "African American"]
    
    # Medication Recommendation Rationale
    print("Based on the patient's profile, here is the recommended three-drug regimen to maximize hypertension treatment:\n")

    # Medication 1: Amlodipine (Calcium Channel Blocker)
    med1 = "Amlodipine"
    reason1 = "A long-acting dihydropyridine calcium channel blocker. This is a first-line agent for African American patients and is highly effective for lowering blood pressure."
    print(f"1. {med1}")
    print(f"   - Rationale: {reason1}\n")

    # Medication 2: Valsartan (Angiotensin II Receptor Blocker - ARB)
    med2 = "Valsartan"
    reason2 = "An ARB is recommended due to the compelling indication of diabetes, as it provides renal protection. ARBs are often preferred over ACE inhibitors in African American patients due to better tolerability."
    print(f"2. {med2}")
    print(f"   - Rationale: {reason2}\n")

    # Medication 3: Spironolactone (Mineralocorticoid Receptor Antagonist - MRA)
    med3 = "Spironolactone"
    reason3 = ("A potassium-sparing diuretic and MRA. This is a key guideline-recommended add-on therapy for resistant hypertension. "
               "It is an effective choice since other first-line diuretics like hydrochlorothiazide and furosemide are not an option for this patient.")
    print(f"3. {med3}")
    print(f"   - Rationale: {reason3}\n")

# Execute the recommendation function
recommend_hypertension_medications()
