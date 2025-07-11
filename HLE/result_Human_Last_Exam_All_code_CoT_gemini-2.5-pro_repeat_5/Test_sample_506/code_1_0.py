def recommend_htn_medications_for_mm():
    """
    This function analyzes the patient's case and recommends a three-drug regimen
    to maximize her hypertension (HTN) treatment.
    """

    # --- Patient Profile Summary ---
    # Age: 27-year-old African American female
    # Condition: Resistant Stage II Hypertension (BP: 145-153/85-89)
    # Comorbidities: Hypercholesterolemia, A1C 6.7%, high Heart Rate (91 bpm)
    # Excluded Medications: Most diuretics, some MRAs, verapamil

    # --- Medication Selection Rationale ---
    # A multi-drug regimen targeting different pathways is required for resistant HTN.

    # Medication 1: Dihydropyridine Calcium Channel Blocker (CCB)
    # Rationale: First-line and highly effective for African American patients per ACC/AHA guidelines.
    medication_1 = "Amlodipine"

    # Medication 2: Angiotensin II Receptor Blocker (ARB)
    # Rationale: First-line agent that inhibits the RAAS system. Beneficial for patients with a high A1C.
    # ARBs are often preferred over ACE inhibitors in African American patients due to lower risk of angioedema.
    medication_2 = "Valsartan"

    # Medication 3: Beta-Blocker with Alpha-Blocking Activity
    # Rationale: An effective add-on therapy for resistant HTN, especially when diuretics are contraindicated.
    # This choice will also help control the patient's elevated heart rate.
    medication_3 = "Carvedilol"

    print("Recommended 3-Medication Regimen for Hypertension:")
    print("1. " + medication_1)
    print("2. " + medication_2)
    print("3. " + medication_3)
    print("\nNote: This regimen combines a CCB, an ARB, and a Beta-Blocker to target multiple pathways, which is suitable for managing resistant hypertension under the given constraints.")

# Run the function to display the recommendations
recommend_htn_medications_for_mm()