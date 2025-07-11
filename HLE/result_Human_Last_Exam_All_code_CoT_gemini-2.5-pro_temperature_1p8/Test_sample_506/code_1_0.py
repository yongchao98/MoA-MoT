def recommend_htn_medications():
    """
    Recommends a 3-drug regimen for the patient's resistant hypertension
    based on her clinical profile and medication restrictions.
    """
    
    # Patient Profile Summary
    patient_info = {
        "Age": 27,
        "Race": "African American",
        "Conditions": ["Resistant Stage II Hypertension", "Diabetes (A1C 6.7%)", "Hypercholesterolemia", "Hypothyroidism"],
        "Blood Pressure": "145-153/85-89 mmHg",
        "BP Goal": "<130/80 mmHg"
    }

    # Medication Rationale
    print("Rationale for Hypertension Regimen:")
    print("-----------------------------------")
    print("The patient has resistant Stage II Hypertension with comorbid diabetes. The goal is to create a potent three-drug combination that follows established guidelines while avoiding her listed medication intolerances.")
    print("\n")
    
    # Recommendation 1: Calcium Channel Blocker
    print("Recommendation 1: Amlodipine (Calcium Channel Blocker)")
    print("   - Rationale: Calcium Channel Blockers are a recommended first-line treatment for hypertension in African American patients. Amlodipine is a dihydropyridine CCB that effectively lowers blood pressure through vasodilation.")
    print("\n")

    # Recommendation 2: ACE Inhibitor
    print("Recommendation 2: Lisinopril (ACE Inhibitor)")
    print("   - Rationale: For patients with both hypertension and diabetes, an ACE inhibitor is recommended to protect the kidneys and control blood pressure. This drug works by inhibiting the Renin-Angiotensin System (RAAS) and combines well with a CCB.")
    print("\n")
    
    # Recommendation 3: Mineralocorticoid Receptor Antagonist
    print("Recommendation 3: Spironolactone (Mineralocorticoid Receptor Antagonist)")
    print("   - Rationale: The patient has resistant hypertension. After a CCB and ACE inhibitor, Spironolactone is the preferred add-on therapy. It works as a potassium-sparing diuretic by blocking aldosterone, which is a key factor in resistant HTN. It is not on the patient's exclusion list and is a critically effective choice here since other diuretics (thiazides and loops) are contraindicated.")
    print("\n")

    print("Final Recommended Regimen:")
    print("1. Amlodipine")
    print("2. Lisinopril")
    print("3. Spironolactone")

# Execute the function to print the recommendations
recommend_htn_medications()