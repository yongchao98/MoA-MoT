def recommend_htn_medications():
    """
    This function outlines the recommended 3-drug regimen for the patient's
    resistant hypertension, considering her comorbidities and medication restrictions.
    """

    # There are no calculations, this is a knowledge-based recommendation.
    # The list below contains the recommended medications and the rationale for each.
    recommendations = {
        "Amlodipine": "This is a dihydropyridine calcium channel blocker (CCB). It is a first-line agent for African American patients and a potent medication for lowering blood pressure.",
        "Losartan": "This is an angiotensin II receptor blocker (ARB). It is a critical component of therapy due to the patient's concurrent type 2 diabetes, as it provides kidney protection.",
        "Chlorthalidone": "This is a thiazide-like diuretic. It is the preferred diuretic in guidelines due to its long duration of action and proven ability to reduce cardiovascular events. It is not on the patient's restricted list and is essential for treating resistant hypertension."
    }

    print("Recommended 3-Medication Regimen to Maximize Hypertension Treatment:")
    print("="*70)
    
    count = 1
    for med, rationale in recommendations.items():
        print(f"{count}. Medication: {med}")
        print(f"   Rationale: {rationale}\n")
        count += 1
        
recommend_htn_medications()
