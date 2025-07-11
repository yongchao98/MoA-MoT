import pandas as pd

def recommend_htn_medication():
    """
    This script analyzes a patient case with resistant hypertension and recommends
    three medications to maximize treatment based on clinical guidelines and a
    list of excluded drugs.
    """
    
    # Patient Profile
    diagnosis = "Resistant Stage II Hypertension"
    current_bp = "145-153/85-89 mmHg"
    goal_bp = "<130/90 mmHg"
    
    # List of medications the patient is unable to take
    excluded_meds = [
        "eplerenone", "finerenone", "hydrochlorothiazide", "indapamide",
        "bumetanide", "furosemide", "torsemide", "metolazone", "verapamil"
    ]
    
    print("Patient Diagnosis:", diagnosis)
    print("Goal BP:", goal_bp)
    print("Excluded Medications:", ', '.join(excluded_meds))
    print("-" * 30)
    print("Recommendation Rationale:\n")

    # --- Medication Recommendations ---
    
    # Medication 1: Chlorthalidone
    rationale_1 = ("1. Rationale for Chlorthalidone:\n"
                   "   - The patient's hypertension is resistant, making a diuretic essential.\n"
                   "   - Common thiazide/thiazide-like diuretics (hydrochlorothiazide, indapamide) are excluded.\n"
                   "   - Chlorthalidone is a potent, long-acting thiazide-like diuretic with strong evidence for "
                   "improving cardiovascular outcomes and is not on the exclusion list.\n")
    rec_1 = "Chlorthalidone"

    # Medication 2: Spironolactone
    rationale_2 = ("2. Rationale for Spironolactone:\n"
                   "   - For resistant hypertension (uncontrolled BP on 3+ drugs), a mineralocorticoid receptor antagonist (MRA) is the recommended fourth-line agent.\n"
                   "   - The MRAs eplerenone and finerenone are excluded, making Spironolactone the appropriate choice.\n"
                   "   - It is highly effective in treating resistant hypertension.\n")
    rec_2 = "Spironolactone"

    # Medication 3: Amlodipine
    rationale_3 = ("3. Rationale for Amlodipine:\n"
                   "   - A dihydropyridine calcium channel blocker (DHP-CCB) is a first-line agent for hypertension.\n"
                   "   - The non-DHP CCB verapamil is excluded, making a DHP-CCB essential.\n"
                   "   - Amlodipine is a potent, long-acting DHP-CCB that effectively lowers BP through vasodilation.\n")
    rec_3 = "Amlodipine"
    
    # Storing recommendations
    recommendations = {
        'Medication': [rec_1, rec_2, rec_3],
        'Rationale': [rationale_1, rationale_2, rationale_3]
    }

    # Print rationale and final list
    print(rationale_1)
    print(rationale_2)
    print(rationale_3)
    
    print("-" * 30)
    print("Final Medication Recommendations to Maximize HTN Treatment:")
    print("-" * 30)
    for med in recommendations['Medication']:
        print(f"- {med}")

# Execute the recommendation function
if __name__ == "__main__":
    recommend_htn_medication()
