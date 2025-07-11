def diagnose_patient():
    """
    Analyzes a clinical case to determine the likely molecular abnormality.
    """
    # Step 1: Define the patient's key clinical data
    patient_age = 18
    patient_percentile_height = "10th - 15th"
    patient_karyotype = "46,XX"
    
    symptoms = [
        "Amenorrhea and infertility",
        "Ovarian dysgenesis (observed on ultrasound)",
        "Persistent short stature ({}th percentile)".format(patient_percentile_height),
        "Fatigue and shortness of breath during physical activity",
        "Intermittent elevated Blood Pressure"
    ]
    
    lab_results = {
        "Karyotype": "{} (Normal chromosomal complement)".format(patient_karyotype)
    }

    # Print the summary of the case
    print("--- Patient Case Summary ---")
    print(f"An {patient_age}-year-old female patient presents with the following:")
    for symptom in symptoms:
        print(f"- {symptom}")
    print("\nKey Lab Result:")
    for test, result in lab_results.items():
        print(f"- {test}: {result}")
        
    # Step 2: Outline the diagnostic reasoning
    print("\n--- Diagnostic Reasoning ---")
    print("1. The symptoms affect multiple systems with high energy demands: ovaries (for fertility), skeletal muscles (for exercise), and the cardiovascular system.")
    print("2. While the phenotype shares features with Turner Syndrome (short stature, ovarian issues), this is ruled out by the normal {} karyotype.".format(patient_karyotype))
    print("3. A problem in cellular energy production would explain the multi-systemic nature of the disease.")
    print("4. The cell's 'powerhouses' are the mitochondria, which contain their own DNA (mtDNA).")
    print("5. A standard karyotype only analyzes the 46 nuclear chromosomes and does not assess mitochondrial DNA.")
    
    # Step 3: State the conclusion
    print("\n--- Conclusion ---")
    conclusion = "A mutation in the mitochondrial DNA (mtDNA)"
    print("Given the involvement of high-energy organ systems and a normal nuclear karyotype, the most probable molecular abnormality is:")
    print(conclusion)

# Execute the diagnostic function
diagnose_patient()