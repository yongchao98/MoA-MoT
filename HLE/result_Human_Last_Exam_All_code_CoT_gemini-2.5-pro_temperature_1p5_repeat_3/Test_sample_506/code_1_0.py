def recommend_htn_medications():
    """
    Recommends a 3-drug regimen for the patient's resistant hypertension based on their clinical profile and medication restrictions.
    """
    
    patient_name = "MM"
    
    # Rationale-driven selection of a 3-drug regimen
    medication_1 = {
        "name": "Amlodipine",
        "class": "Dihydropyridine Calcium Channel Blocker (CCB)",
        "rationale": "First-line agent for hypertension, particularly effective in African American patients. It is not on the patient's exclusion list."
    }
    
    medication_2 = {
        "name": "Losartan",
        "class": "Angiotensin II Receptor Blocker (ARB)",
        "rationale": "Recommended for patients with diabetes (A1C 6.7%) to provide renal protection. Preferred over an ACE inhibitor due to a lower risk of cough."
    }
    
    medication_3 = {
        "name": "Chlorthalidone",
        "class": "Thiazide-like Diuretic",
        "rationale": "A potent, long-acting diuretic highly effective for resistant hypertension. It is not on the patient's exclusion list and is often preferred over hydrochlorothiazide."
    }
    
    recommendations = [medication_1, medication_2, medication_3]
    
    print(f"Based on the clinical profile for patient {patient_name}, here are 3 recommended medications to maximize hypertension treatment:\n")
    
    for i, med in enumerate(recommendations, 1):
        print(f"Recommendation {i}:")
        print(f"  Medication: {med['name']}")
        print(f"  Drug Class: {med['class']}")
        print(f"  Rationale: {med['rationale']}\n")

if __name__ == "__main__":
    recommend_htn_medications()