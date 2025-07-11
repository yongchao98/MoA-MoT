def diagnose_patient():
    """
    Analyzes the patient's clinical data to arrive at a diagnosis.
    """
    # Step 1: Define patient's key data points from the case
    patient_age = 64
    patient_bmi = 39
    is_smoker = True
    years_smoking = 15
    lesion_locations = ["axillary folds", "inframammary folds", "inguinal regions"]
    lesion_types = ["large bullae", "erythematous plaques", "purulent nodules"]

    # Step 2: Print the analysis and reasoning
    print("--- Diagnostic Analysis ---")
    print(f"Patient Profile: A {patient_age}-year-old woman with a BMI of {patient_bmi}.")
    print(f"Significant Risk Factors: Obesity (BMI > 30) and a {years_smoking}-year smoking history.")
    print("\nClinical Findings:")
    print(f"Locations: Lesions are present in multiple skin folds: {', '.join(lesion_locations)}.")
    print(f"Morphology: The key finding is the presence of 'purulent nodules', especially in the inguinal regions.")

    print("\nReasoning:")
    print("1. The combination of lesion locations (axillary, inguinal) is characteristic of a disease process affecting apocrine glands.")
    print("2. The presence of 'purulent nodules' is a hallmark sign of Hidradenitis Suppurativa (HS).")
    print("3. The patient's major risk factors, obesity (BMI 39) and smoking, are strongly associated with the development and severity of HS.")
    print("4. While other conditions can affect skin folds, none of the other choices typically present with deep, inflammatory, purulent nodules in this distribution.")

    print("\n--- Conclusion ---")
    print("The clinical picture of purulent nodules in the axillary and inguinal folds in an obese patient who smokes is classic for Hidradenitis Suppurativa.")

# Execute the diagnostic function
diagnose_patient()