def recommend_htn_medications():
    """
    Analyzes the patient's case and recommends a 3-drug regimen for resistant hypertension.
    
    Patient Profile:
    - 27-year-old African American female
    - Resistant Stage II Hypertension (BP: 145-153/85-89)
    - Comorbidities: Hypercholesterolemia, Diabetes (A1C 6.7%), undertreated Hypothyroidism
    - Excluded Medications: eplerenone, finerenone, hydrochlorothiazide, indapamide, 
      bumetanide, furosemide, torsemide, metolazone, and verapamil.
    """
    
    # Medication recommendations based on clinical guidelines for resistant hypertension
    # in an African American patient, considering the exclusion list.
    
    medication_1 = {
        "name": "Amlodipine",
        "class": "Dihydropyridine Calcium Channel Blocker (CCB)",
        "rationale": "A first-line agent for hypertension, particularly effective in African American patients. It is not on the exclusion list."
    }
    
    medication_2 = {
        "name": "Olmesartan",
        "class": "Angiotensin II Receptor Blocker (ARB)",
        "rationale": "Blocks the renin-angiotensin system. ARBs are preferred for this patient to avoid the cough associated with ACE inhibitors and are beneficial in patients with diabetes."
    }
    
    medication_3 = {
        "name": "Spironolactone",
        "class": "Mineralocorticoid Receptor Antagonist (MRA) / Potassium-Sparing Diuretic",
        "rationale": "The recommended add-on therapy for resistant hypertension. It is a crucial choice as other diuretics and MRAs are on the patient's exclusion list."
    }
    
    recommendations = [medication_1, medication_2, medication_3]
    
    print("Based on the patient's profile and medication constraints, here are three recommended medications to maximize hypertension treatment:\n")
    
    for i, med in enumerate(recommendations, 1):
        print(f"Recommendation {i}:")
        print(f"  - Medication: {med['name']}")
        print(f"  - Rationale: {med['rationale']}\n")

# Execute the function to print the recommendations
recommend_htn_medications()