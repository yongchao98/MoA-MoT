def recommend_htn_medications():
    """
    This function recommends three medications to maximize hypertension (HTN) treatment
    for the patient MM based on her clinical profile and medication restrictions.
    """
    # Patient Profile: 27-year-old African American female with resistant Stage II HTN,
    # Type 2 Diabetes (A1C=6.7%), and hypercholesterolemia.
    # Excluded medications prevent the use of most standard diuretics.

    # Medication 1: Amlodipine (Calcium Channel Blocker)
    # Rationale: First-line recommendation for African American patients. It is a potent
    # dihydropyridine CCB, and this class is not on the patient's exclusion list.
    medication_1 = "Amlodipine"

    # Medication 2: Losartan (Angiotensin II Receptor Blocker - ARB)
    # Rationale: Strong recommendation for patients with diabetes to provide renal
    # protection. It effectively lowers blood pressure by blocking the RAAS system.
    medication_2 = "Losartan"

    # Medication 3: Spironolactone (Mineralocorticoid Receptor Antagonist - MRA)
    # Rationale: Recommended add-on therapy for resistant hypertension. It acts as a
    # potassium-sparing diuretic. It is notably NOT on the patient's exclusion list,
    # making it the ideal third agent where other diuretics are contraindicated.
    medication_3 = "Spironolactone"
    
    print("Based on the patient's profile and medication restrictions, the following 3 medications are recommended to maximize her hypertension treatment:")
    print(f"1. {medication_1}")
    print(f"2. {medication_2}")
    print(f"3. {medication_3}")
    
    print("\nImportant Note: This combination requires monitoring of serum potassium levels due to the risk of hyperkalemia from Losartan and Spironolactone.")

recommend_htn_medications()