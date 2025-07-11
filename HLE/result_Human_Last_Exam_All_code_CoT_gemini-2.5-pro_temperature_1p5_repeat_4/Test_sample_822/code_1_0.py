import sys

def diagnose_case_from_vignette():
    """
    Analyzes a clinical vignette to determine the most likely underlying disease.
    """
    
    # Define key findings from the clinical case
    patient_profile = {
        "Age": "62-year-old male",
        "Risk Factors": ["20-pack-year smoking", "Ship building (asbestos exposure risk)"]
    }

    chronic_symptoms = {
        "Description": "Initial multi-system presentation",
        "Symptoms": ["Fatigue", "Joint pain/swelling (wrists, ankles, elbows)", "Neurological (dizziness, confusion)", "Systemic (bruising, difficulty swallowing, loss of appetite)"],
        "Imaging": "Multiple pulmonary nodules on Chest X-ray"
    }

    treatment_and_complication = {
        "Initial Treatment": "Steroids and NSAIDs",
        "Consequence": "Steroid use leads to an immunocompromised state.",
        "Acute Event": "Development of fever, productive cough, shortness of breath, cutaneous lesions.",
        "Failed Treatment": "Aminoglycoside therapy was ineffective.",
        "Outcome": "Death from septic shock"
    }

    # Step-by-step diagnostic reasoning
    print("Step 1: Analyzing Patient History and Chronic Symptoms")
    print(f"The patient's profile and initial symptoms ({chronic_symptoms['Symptoms']} plus {chronic_symptoms['Imaging']}) strongly suggest a systemic autoimmune vasculitis.")
    print("-" * 30)
    
    print("Step 2: Forming a Differential Diagnosis")
    print("The combination of upper/lower respiratory tract involvement (nodules, shortness of breath) and systemic features (joint pain, fatigue) makes Granulomatosis with Polyangiitis (GPA) a primary diagnosis.")
    print("-" * 30)

    print("Step 3: Evaluating the Terminal Illness")
    print(f"The patient was started on steroids, which caused immunosuppression. He then developed a severe infection.")
    print(f"The key clue is that '{treatment_and_complication['Failed Treatment']}'. This points towards an atypical opportunistic pathogen.")
    print("Nocardiosis is a classic opportunistic infection that causes lung and skin lesions in patients with GPA on steroid therapy, and it does not respond to aminoglycosides.")
    print("-" * 30)

    print("Step 4: Synthesizing a Final Diagnosis")
    print("The entire clinical picture, from the initial presentation to the fatal complication, is best explained by a primary underlying disease.")
    final_diagnosis = "Granulomatosis with Polyangiitis (GPA)"
    print(f"The most likely underlying disease the patient experienced was: {final_diagnosis}")
    
    # The final answer format as requested by the user
    # Suppressing this part from the standard output to only have it at the very end.
    original_stdout = sys.stdout
    sys.stdout = open('/dev/null', 'w')
    print(f'<<<{final_diagnosis}>>>')
    sys.stdout.close()
    sys.stdout = original_stdout


diagnose_case_from_vignette()
<<<Granulomatosis with Polyangiitis (GPA)>>>