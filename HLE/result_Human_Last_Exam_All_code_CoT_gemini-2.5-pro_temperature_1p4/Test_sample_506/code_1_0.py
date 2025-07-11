def recommend_hypertension_medications():
    """
    This script recommends a 3-drug regimen for a patient with resistant hypertension,
    diabetes, and other comorbidities, based on a clinical case study.
    """

    print("### Rationale for Recommended Hypertension Regimen ###\n")
    print("The patient has resistant Stage II Hypertension and comorbid Diabetes. The goal is to recommend an optimal 3-drug regimen to maximize blood pressure control, considering her specific profile and medication constraints.\n")
    
    # Medication 1: Calcium Channel Blocker (CCB)
    medication_1 = "Amlodipine"
    rationale_1 = ("A long-acting dihydropyridine Calcium Channel Blocker (CCB). "
                   "CCBs are a first-line therapy for hypertension and are highly effective in African American patients. "
                   "It is not on the patient's exclusion list.")
    
    # Medication 2: Angiotensin II Receptor Blocker (ARB)
    medication_2 = "Losartan"
    rationale_2 = ("An Angiotensin II Receptor Blocker (ARB). Given the patient's diagnosis of diabetes (A1C 6.7%), "
                   "an ARB or ACE inhibitor is essential for renal protection. ARBs are often chosen for their favorable side-effect profile.")
                   
    # Medication 3: Thiazide-like Diuretic
    medication_3 = "Chlorthalidone"
    rationale_3 = ("A long-acting thiazide-like diuretic. A diuretic is a cornerstone of multi-drug hypertension therapy. "
                   "While many common diuretics are on the exclusion list, Chlorthalidone is not. It is more potent and longer-acting "
                   "than HCTZ, with strong data showing improved cardiovascular outcomes.")

    recommendations = {
        medication_1: rationale_1,
        medication_2: rationale_2,
        medication_3: rationale_3
    }
    
    print("Based on this analysis, the following 3 medications are recommended:\n")
    
    count = 1
    for med, rationale in recommendations.items():
        print(f"{count}. Medication: {med}")
        print(f"   Rationale: {rationale}\n")
        count += 1
        
    print("-" * 60)
    print("This combination creates a powerful, synergistic effect by targeting volume, vascular resistance, and the renin-angiotensin system, and is aligned with current clinical guidelines.")

if __name__ == "__main__":
    recommend_hypertension_medications()