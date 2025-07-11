def solve_clinical_case():
    """
    This function analyzes the provided clinical case and prints the likely diagnosis.
    """

    patient_age = 62
    smoking_history_pack_years = 20
    
    # Key clinical findings
    initial_symptoms = ["polyarthritis", "fatigue", "multiple pulmonary nodules"]
    progression_symptoms = ["confusion", "cutaneous lesions", "productive cough"]
    immunosuppression = "Started on steroids"
    failed_treatment = "Aminoglycoside therapy"
    
    # Reasoning
    print("Clinical Reasoning:")
    print(f"1. An {patient_age}-year-old man with risk factors (a {smoking_history_pack_years}-pack-year history) was treated with steroids, leading to immunosuppression.")
    print(f"2. He developed a severe disseminated infection with symptoms affecting the lungs ({initial_symptoms[2]}), skin ({progression_symptoms[1]}), and likely the central nervous system ({progression_symptoms[0]}).")
    print(f"3. The infection proved resistant to '{failed_treatment}', which is a major clue.")
    print("4. Nocardia is a classic opportunistic pathogen in steroid-treated patients and is characteristically not susceptible to aminoglycosides.")
    
    diagnosis = "Nocardiosis"
    
    print("\nConclusion:")
    print(f"The combination of the patient's immunocompromised state, multi-systemic infection (pulmonary, cutaneous, CNS), and failure to respond to aminoglycoside therapy strongly indicates the disease was {diagnosis}.")

solve_clinical_case()