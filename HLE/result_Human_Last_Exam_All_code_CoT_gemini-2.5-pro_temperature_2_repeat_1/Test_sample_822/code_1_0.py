def diagnose_patient():
    """
    Analyzes clinical findings to suggest the most likely disease.
    This function uses a simplified scoring model to illustrate the diagnostic logic.
    """
    # Key clinical findings from the case
    findings = {
        "multi_system_inflammation": True, # Arthritis, fatigue, bruising
        "multiple_pulmonary_nodules": True,
        "opportunistic_infection_pattern": True, # Pulmonary + Cutaneous, resistant to aminoglycosides
        "high_cancer_risk_factors": True # Smoking, asbestos exposure
    }

    # Potential diagnoses and how well they fit the findings
    disease_scores = {
        "Granulomatosis with Polyangiitis (GPA)": 0,
        "Lung Cancer w/ Paraneoplastic Syndrome": 0,
        "Sarcoidosis": 0,
        "Tuberculosis": 0
    }

    # Scoring logic
    # GPA
    if findings["multi_system_inflammation"]:
        disease_scores["Granulomatosis with Polyangiitis (GPA)"] += 1
    if findings["multiple_pulmonary_nodules"]:
        disease_scores["Granulomatosis with Polyangiitis (GPA)"] += 1
    if findings["opportunistic_infection_pattern"]:
        # Nocardia is a classic opportunistic infection in GPA
        disease_scores["Granulomatosis with Polyangiitis (GPA)"] += 1
    
    # Lung Cancer
    if findings["high_cancer_risk_factors"]:
        disease_scores["Lung Cancer w/ Paraneoplastic Syndrome"] += 1
    if findings["multiple_pulmonary_nodules"]:
        # Could be metastases
        disease_scores["Lung Cancer w/ Paraneoplastic Syndrome"] += 1
    if findings["multi_system_inflammation"]:
         # Paraneoplastic syndromes can mimic this
        disease_scores["Lung Cancer w/ Paraneoplastic Syndrome"] += 1

    # Sarcoidosis
    if findings["multi_system_inflammation"]:
        disease_scores["Sarcoidosis"] += 1
    if findings["multiple_pulmonary_nodules"]:
        disease_scores["Sarcoidosis"] += 1
        
    # Tuberculosis
    if findings["multiple_pulmonary_nodules"]:
        disease_scores["Tuberculosis"] += 1
        
    print("Clinical Feature Matching Score:")
    for disease, score in disease_scores.items():
        print(f"- {disease}: {score}")

    # Determine the most likely diagnosis
    most_likely_disease = max(disease_scores, key=disease_scores.get)

    print("\nConclusion:")
    print("The constellation of symptoms, including multi-system inflammation (arthritis),")
    print("multiple pulmonary nodules, and a subsequent fatal opportunistic infection")
    print("with a specific pattern (pulmonary and cutaneous lesions, resistant to aminoglycosides)")
    print("is most characteristic of Granulomatosis with Polyangiitis (GPA).")
    
    # Returning the final answer in the specified format
    print(f"\nFinal Answer: {most_likely_disease}")

diagnose_patient()