def recommend_htn_medications():
    """
    This function recommends three medications for a patient with resistant hypertension,
    based on their clinical profile and medication restrictions.
    """

    # Patient Profile: 27-year-old African American female with resistant HTN, T2DM.
    # Restrictions: eplerenone, finerenone, HCTZ, indapamide, loop diuretics, verapamil.

    # 1. Amlodipine (CCB): First-line for African American patients.
    medication1 = "Amlodipine"

    # 2. Losartan (ARB): For renal protection due to T2DM (compelling indication).
    medication2 = "Losartan"

    # 3. Spironolactone (MRA): Guideline-recommended for resistant HTN, especially
    # when other diuretics are not an option. It is not on the restricted list.
    medication3 = "Spironolactone"

    recommendations = [medication1, medication2, medication3]

    print("Recommended 3-drug regimen for maximizing hypertension treatment:")
    for med in recommendations:
        print(f"- {med}")

recommend_htn_medications()