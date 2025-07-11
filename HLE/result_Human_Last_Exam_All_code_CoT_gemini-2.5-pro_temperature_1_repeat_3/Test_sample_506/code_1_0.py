def recommend_hypertension_medications():
    """
    Analyzes the patient's case and recommends a 3-drug regimen
    to maximize hypertension treatment.
    """
    # Based on the analysis, the three recommended medications are:
    # 1. Losartan (ARB) - For RAAS inhibition and diabetic kidney protection.
    # 2. Amlodipine (CCB) - Highly effective in African American patients.
    # 3. Spironolactone (MRA) - Guideline-recommended for resistant hypertension.
    
    recommended_medications = ["Losartan", "Amlodipine", "Spironolactone"]
    
    print("Based on the patient's profile of resistant hypertension and co-morbid conditions, the following 3 medications are recommended to form an optimal treatment regimen:")
    print("-" * 70)
    
    # Print each medication with its rationale
    print(f"1. Medication: {recommended_medications[0]}")
    print("   Rationale: This is an Angiotensin Receptor Blocker (ARB). It is recommended for patients with diabetes to protect the kidneys and is highly effective for hypertension. It is chosen over an ACE-inhibitor to reduce the risk of cough.\n")
    
    print(f"2. Medication: {recommended_medications[1]}")
    print("   Rationale: This is a Dihydropyridine Calcium Channel Blocker (CCB), which is a first-line and highly effective treatment for hypertension in African American patients.\n")
    
    print(f"3. Medication: {recommended_medications[2]}")
    print("   Rationale: This is a Mineralocorticoid Receptor Antagonist (MRA) and a key medication for resistant hypertension. It addresses blood pressure elevation caused by the hormone aldosterone and is not on the patient's restricted list.\n")
          
    print("-" * 70)
    print("This combination targets three different mechanisms to provide maximum blood pressure control.")
    print("Monitoring of blood pressure, potassium levels, and kidney function is essential after starting this regimen.")

# Execute the function to print the recommendation
recommend_hypertension_medications()