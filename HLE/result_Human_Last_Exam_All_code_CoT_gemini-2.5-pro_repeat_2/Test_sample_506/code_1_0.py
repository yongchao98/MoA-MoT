def recommend_htn_medications():
    """
    Recommends a three-drug regimen for the patient's resistant hypertension based on
    current guidelines and patient-specific factors.
    """
    print("Based on the patient's profile (African American female with resistant Stage II HTN) and specified medication intolerances, here is a recommended 3-drug regimen to maximize blood pressure control:")
    print("-" * 80)

    # Recommendation 1: Dihydropyridine Calcium Channel Blocker (CCB)
    print("1. Amlodipine")
    print("   Rationale: A long-acting dihydropyridine CCB is a first-line, preferred agent for hypertension in African American patients. It is a potent vasodilator and is not on the patient's exclusion list.")
    print("-" * 80)

    # Recommendation 2: Angiotensin II Receptor Blocker (ARB)
    print("2. Losartan")
    print("   Rationale: An ARB blocks the renin-angiotensin system and is a standard part of combination therapy. This class is also beneficial given the patient's A1C of 6.7%. An ARB is chosen over an ACE-Inhibitor to reduce the risk of cough.")
    print("-" * 80)

    # Recommendation 3: Mineralocorticoid Receptor Antagonist (MRA)
    print("3. Spironolactone")
    print("   Rationale: For resistant hypertension, especially when standard diuretics are not an option, a mineralocorticoid receptor antagonist is the guideline-recommended add-on therapy. Spironolactone targets aldosterone, a key driver of resistant HTN, and is not on the patient's exclusion list.")
    print("-" * 80)

if __name__ == '__main__':
    recommend_htn_medications()