def diagnose_patient():
    """
    This function analyzes the patient's clinical data and determines the most likely diagnosis.
    """
    # Patient Data
    patient_age = 57
    symptoms = ["dyspnea", "chronic cough", "acid reflux"]
    history = ["COPD"]
    imaging_findings = {"CT Chest": "mass of the vertebrae"}
    lab_findings = {"Creatinine": 2.1} # Elevated, indicates renal dysfunction

    print("Analyzing the clinical case based on the provided information...")
    print("-" * 50)
    print(f"Patient is a {patient_age}-year-old woman with a history of {history[0]}.")
    print(f"She presents with symptoms: {', '.join(symptoms)}.")
    print(f"A key finding from imaging is a {imaging_findings['CT Chest']}.")
    print("-" * 50)
    print("Evaluating the potential diagnoses:\n")

    # Dictionary to hold the reasoning for each choice
    analysis = {
        "A. Aspiration pneumonitis": "Incorrect. Does not explain the chronic symptoms or the vertebral mass.",
        "B. Aspiration pneumonia": "Incorrect. This is an infection and does not explain the vertebral mass.",
        "C. Achalasia": "Incorrect. This is an esophageal motility disorder and does not cause a bone mass.",
        "D. Adenocarcinoma": "Correct. This is a type of cancer. A primary lung adenocarcinoma explains the respiratory symptoms (cough, dyspnea) and also commonly metastasizes to bone, which explains the vertebral mass. This provides a single unifying diagnosis for the most critical findings.",
        "E. COPD": "Incorrect. This is a pre-existing condition and does not cause a new mass in the vertebrae."
    }

    for option, reason in analysis.items():
        print(f"Diagnosis Option: {option}")
        print(f"Reasoning: {reason}\n")

    print("-" * 50)
    print("Conclusion: The most compelling evidence is the vertebral mass, which strongly suggests metastatic cancer. Adenocarcinoma is the only option that accounts for all the patient's major findings.")
    print("Final Answer: D. Adenocarcinoma")

# Execute the diagnostic analysis
diagnose_patient()