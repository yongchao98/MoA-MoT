import sys
# The user wants me to act as an AI assistant to recommend 3 medications for resistant hypertension.
# I need to analyze the patient's case, consider her comorbidities and the list of medications she cannot take.
# Finally, I will print out the reasoning and the recommended medications in the specified format.

def recommend_htn_medications():
    """
    Analyzes the patient case and recommends a 3-drug regimen for resistant hypertension.
    """
    # Patient Profile Summary
    patient_profile = {
        "Age": 27,
        "Ethnicity": "African American",
        "Condition": "Resistant Stage II Hypertension",
        "BP": "145-153/85-89 mmHg",
        "Goal BP": "<130/90 mmHg",
        "Comorbidities": ["Type 2 Diabetes (A1C 6.7%)", "Hypercholesterolemia", "Undertreated Hypothyroidism"],
        "Excluded Medications": ["eplerenone", "finerenone", "hydrochlorothiazide", "indapamide", 
                               "bumetanide", "furosemide", "torsemide", "metolazone", "verapamil"]
    }

    # Recommendation Rationale
    print("Based on the patient's profile of resistant hypertension, Type 2 Diabetes, and African American ethnicity, here is a recommended 3-drug regimen to maximize treatment while respecting her medication restrictions:\n")

    # Medication 1
    print("1. Losartan (an Angiotensin II Receptor Blocker - ARB)")
    print("   - Rationale: First-line therapy for patients with both hypertension and diabetes due to its proven renal-protective benefits. ARBs are chosen over ACE-Inhibitors to avoid the potential side effect of a persistent cough.\n")

    # Medication 2
    print("2. Amlodipine (a Dihydropyridine Calcium Channel Blocker - CCB)")
    print("   - Rationale: CCBs are highly effective in African American patients and are a first-line agent. Combining an ARB with a CCB is a powerful and synergistic strategy for blood pressure control. This class of CCB is permitted as only verapamil is excluded.\n")

    # Medication 3
    print("3. Spironolactone (a Mineralocorticoid Receptor Antagonist - MRA)")
    print("   - Rationale: The patient has resistant hypertension, and all standard diuretics are excluded. Spironolactone is the guideline-recommended add-on for resistant HTN. It functions as a potassium-sparing diuretic and targets the aldosterone pathway, which is often implicated in difficult-to-treat hypertension. While other MRAs are excluded, Spironolactone is a viable and potent option.\n")
    
    # Final answer format
    final_answer = "<<<Losartan, Amlodipine, Spironolactone>>>"
    # The following line should not be commented out in the final output.
    # It prints the final answer as requested by the user.
    # Using sys.stdout.write to avoid adding an extra newline if print does.
    sys.stdout.write(final_answer)

recommend_htn_medications()