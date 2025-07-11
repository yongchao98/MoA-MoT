def diagnose_patient():
    """
    Analyzes the patient's case to determine the most likely diagnosis.
    """
    # Patient Data
    age = 62  # years
    smoking_history_pack_years = 20

    # Key Clinical Findings
    risk_factors = [
        "Age: {} years".format(age),
        "Smoking: {} pack-year history".format(smoking_history_pack_years),
        "Occupation: Ship building (asbestos exposure risk)"
    ]

    symptoms_and_signs = [
        "Initial: Polyarthritis (wrists, ankles, elbows), fatigue",
        "Progression: Dizziness, confusion, bruising, dysphagia, weight loss",
        "Imaging: Multiple pulmonary nodules on chest X-ray",
        "Terminal Event: Severe opportunistic infection leading to septic shock"
    ]

    # Diagnostic Reasoning
    print("Patient Profile and Risk Factors:")
    for factor in risk_factors:
        print(f"- {factor}")
    print("\nClinical Presentation Summary:")
    for sign in symptoms_and_signs:
        print(f"- {sign}")

    print("\nConclusion:")
    print("The combination of major risk factors (smoking and asbestos exposure), a clinical picture classic for paraneoplastic syndromes (arthritis, neurologic symptoms), and imaging showing multiple pulmonary nodules points overwhelmingly to a single diagnosis.")
    print("The final infectious episode represents a common complication in an immunocompromised patient with this underlying disease.")
    
    # Final Diagnosis
    final_diagnosis = "Lung Cancer"
    print(f"\nThe most likely diagnosis is: {final_diagnosis}")

diagnose_patient()