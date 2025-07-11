def diagnose_patient():
    """
    This script analyzes the patient's clinical data to determine the most likely diagnosis.
    """

    # Patient Data from the case
    age = 57
    creatinine = 2.1
    symptoms = ["dyspnea", "chronic cough", "acid reflux"]
    history = "COPD"
    imaging_finding = "mass of the vertebrae on CT scan"

    # Printing the clinical information and the numbers mentioned
    print("--- Patient Clinical Summary ---")
    print(f"Age: {age}")
    print(f"Significant Lab Value (Creatinine): {creatinine}")
    print(f"Key Imaging Finding: {imaging_finding}")
    print("\n--- Diagnostic Reasoning ---")
    print("The primary symptoms like dyspnea and chronic cough are consistent with the patient's history of COPD.")
    print("However, the most significant and differentiating finding is the 'mass of the vertebrae'.")
    print("This finding is a strong indicator of metastatic cancer, which is cancer that has spread to the bones.")
    print("Reviewing the options:")
    print(" - Aspiration, Achalasia, and COPD do not cause bone masses.")
    print(" - Adenocarcinoma is a type of cancer (often lung cancer) that explains both the respiratory symptoms and the metastatic bone mass.")
    print("\n--- Final Conclusion ---")
    print("The diagnosis that best explains the combination of respiratory symptoms and a vertebral mass is Adenocarcinoma.")

diagnose_patient()