def solve_case_study():
    """
    Analyzes the patient's case study to determine the most likely disease.
    """
    patient_profile = {
        "Age": 62,
        "Smoking History (pack-years)": 20,
        "Occupation": "Ship building (risk for inhaled dust like silica/asbestos)"
    }

    chronic_symptoms = [
        "Fatigue",
        "Polyarthritis (wrists, ankles, elbows)",
        "Shortness of breath",
        "Systemic symptoms (dizziness, confusion, bruising, dysphagia)",
        "Chest X-ray: Multiple pulmonary nodules"
    ]

    acute_illness_symptoms = [
        "Fever",
        "Productive cough (green sputum)",
        "Cutaneous lesions",
        "Refractory to Aminoglycoside therapy",
        "Outcome: Death from septic shock"
    ]

    disease_analysis = {
        "Key Clue 1": "The combination of chronic respiratory symptoms (fatigue, shortness of breath) and pulmonary nodules points to an underlying lung pathology. Pulmonary Alveolar Proteinosis (PAP) is a strong possibility.",
        "Key Clue 2": "The terminal illness with lung, skin (cutaneous lesions), and likely CNS (confusion) involvement is classic for disseminated Nocardiosis.",
        "Key Clue 3": "PAP creates a protein-rich environment in the lungs, making patients highly susceptible to opportunistic infections, especially Nocardia.",
        "Key Clue 4": "Nocardia infections often do not respond to single-agent therapies like aminoglycosides and require specific, long-term combination antibiotics."
    }

    conclusion = "Pulmonary Alveolar Proteinosis (PAP) complicated by disseminated Nocardiosis"

    print("Patient Profile:")
    for key, value in patient_profile.items():
        print(f"- {key}: {value}")
    
    print("\nInitial Chronic Symptoms:")
    for symptom in chronic_symptoms:
        print(f"- {symptom}")

    print("\nFinal Acute Illness:")
    for finding in acute_illness_symptoms:
        print(f"- {finding}")

    print("\nDiagnostic Reasoning:")
    print(disease_analysis["Key Clue 1"])
    print(disease_analysis["Key Clue 2"])
    print(disease_analysis["Key Clue 3"])
    print(disease_analysis["Key Clue 4"])
    
    print("\n---")
    print(f"Final Conclusion: The most likely disease experienced by the patient was {conclusion}.")

solve_case_study()
<<<Pulmonary Alveolar Proteinosis (PAP) complicated by disseminated Nocardiosis>>>