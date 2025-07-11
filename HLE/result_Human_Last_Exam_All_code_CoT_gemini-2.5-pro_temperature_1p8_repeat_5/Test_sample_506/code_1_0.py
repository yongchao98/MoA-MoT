def recommend_htn_medications():
    """
    Analyzes the patient case and recommends a 3-drug regimen to maximize hypertension treatment.
    """
    
    patient_case_summary = """
    Patient: 27-year-old African American female.
    Key Conditions: Resistant Stage II Hypertension (BP 145-153/85-89), Diabetes (A1C 6.7%), Hypercholesterolemia.
    Goal: Maximize hypertension treatment with three medications, considering a list of medications she cannot take.
    """
    
    rationale = """
    Rationale:
    Based on the patient's profile and current hypertension guidelines, the most effective approach is a combination of three first-line agents with complementary mechanisms of action. This regimen is selected to maximize blood pressure reduction while also providing benefits for her diabetes. The chosen medications are not on the patient's exclusion list.
    """
    
    recommendations = {
        "1. Losartan": "This is an Angiotensin II Receptor Blocker (ARB). It is a highly effective blood pressure medication and is strongly recommended for patients with diabetes to help protect their kidneys. It is generally better tolerated than an ACE inhibitor.",
        "2. Amlodipine": "This is a dihydropyridine Calcium Channel Blocker (CCB). CCBs are a first-line and particularly effective choice for treating hypertension in African American patients. It works by relaxing blood vessels.",
        "3. Chlorthalidone": "This is a thiazide-like diuretic. While many diuretics are on the patient's exclusion list, this one is not. It is a potent, long-acting diuretic that is recommended by guidelines for its superior 24-hour blood pressure control and proven ability to reduce cardiovascular events."
    }
    
    print("Recommended 3-Drug Regimen for Maximal HTN Treatment:")
    print("-" * 55)
    for med, reason in recommendations.items():
        print(f"{med}: {reason}\n")

# Execute the function to print the recommendations
recommend_htn_medications()