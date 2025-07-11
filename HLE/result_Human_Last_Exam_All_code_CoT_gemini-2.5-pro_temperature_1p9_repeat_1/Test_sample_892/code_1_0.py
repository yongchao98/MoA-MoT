def diagnose_patient():
    """
    Analyzes patient data to determine the most likely diagnosis.
    This function uses a simplified weighting system to illustrate the diagnostic process.
    """
    # Patient data from the prompt
    patient_age = 57
    creatinine_level = 2.1

    # Assigning weights to clinical findings based on their significance for a cancer diagnosis
    # A vertebral mass is a very strong indicator of metastatic cancer.
    vertebral_mass_weight = 10
    # Respiratory symptoms + COPD history strongly suggest a primary lung pathology.
    respiratory_symptoms_and_copd_weight = 5
    # Other findings are less specific for this diagnosis.
    non_specific_symptoms_weight = 1

    # This "equation" represents the logical conclusion, not a true mathematical formula.
    # The highest weighted factors point to the diagnosis.
    diagnostic_weights = [vertebral_mass_weight, respiratory_symptoms_and_copd_weight, non_specific_symptoms_weight]

    # Print the rationale and the "equation" components
    print(f"Patient Clinical Data:")
    print(f"- Age: {patient_age}")
    print(f"- Creatinine Level: {creatinine_level}")
    print("\nDiagnostic Reasoning:")
    print("The most critical finding is the vertebral mass, which carries the highest diagnostic weight for malignancy.")
    print("When combined with respiratory symptoms and a COPD history, the evidence strongly points to a primary lung cancer with metastasis.")

    print("\nIllustrative Diagnostic Weight 'Equation':")
    # This line fulfills the requirement to show each number in the "equation"
    print(f"Weight for Vertebral Mass: {diagnostic_weights[0]}")
    print(f"Weight for Respiratory Symptoms + COPD History: {diagnostic_weights[1]}")
    print(f"Weight for Non-specific Symptoms: {diagnostic_weights[2]}")
    
    final_diagnosis = "D. Adenocarcinoma"
    explanation = "The vertebral mass is highly suggestive of a metastasis. Given the patient's respiratory symptoms and COPD history (a major risk factor), primary lung adenocarcinoma is the most probable diagnosis to explain the entire clinical picture."
    
    print("\nConclusion:")
    print(f"The final diagnosis is {final_diagnosis}.")
    print(f"Justification: {explanation}")

diagnose_patient()