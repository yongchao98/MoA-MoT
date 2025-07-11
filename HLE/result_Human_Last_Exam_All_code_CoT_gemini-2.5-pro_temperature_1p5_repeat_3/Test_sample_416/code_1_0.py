def diagnose_patient_condition():
    """
    Analyzes a patient's clinical findings to determine the most likely diagnosis.
    """
    # Key findings from the patient case
    patient_findings = {
        "presentation": "Acute red, swollen, painful ankle",
        "imaging": "X-rays negative for acute abnormality",
        "synovial_fluid_crystals": "None",
        "synovial_fluid_organisms_wbc": "None",
        "treatment_response": "Worsened despite NSAIDs and steroids"
    }

    print("Analyzing patient findings:")
    for key, value in patient_findings.items():
        print(f"- {key.replace('_', ' ').title()}: {value}")

    print("\nDiagnostic Reasoning:")

    # Rule out diagnoses based on definitive tests
    if patient_findings["synovial_fluid_crystals"] == "None":
        print("- Absence of crystals in synovial fluid rules out Gout and Pseudogout.")
    
    if patient_findings["synovial_fluid_organisms_wbc"] == "None":
        print("- Absence of organisms or white blood cells in synovial fluid makes Septic Arthritis highly unlikely.")

    # Consider remaining diagnoses
    print("- Osteoarthritis is less likely due to the highly inflammatory presentation and lack of response to steroids.")
    print("- Charcot Arthropathy is highly suspected due to the combination of an acutely inflamed joint, negative initial X-rays, non-diagnostic fluid, and worsening symptoms despite anti-inflammatory treatment. This presentation is classic for the acute phase of neuropathic destruction.")
    
    final_diagnosis = "Charcot Arthropathy"
    
    print(f"\nConclusion: The clinical picture strongly points towards a diagnosis of {final_diagnosis}.")

# Execute the diagnostic logic
diagnose_patient_condition()