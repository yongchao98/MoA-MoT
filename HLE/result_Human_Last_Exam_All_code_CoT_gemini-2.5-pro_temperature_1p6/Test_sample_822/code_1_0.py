def diagnose_patient():
    """
    Analyzes patient symptoms and history to suggest a likely diagnosis.
    """
    # Potential diagnoses and their scores, initialized to 0
    diagnoses = {
        "Pulmonary Nocardiosis with Dissemination": 0,
        "Granulomatosis with Polyangiitis (GPA)": 0,
        "Metastatic Lung Cancer": 0,
        "Tuberculosis": 0
    }

    # Clinical findings from the case study and their impact on diagnosis scores
    # This simulates a differential diagnosis process.
    findings = {
        "Immunocompromised state (steroid use)": {"Pulmonary Nocardiosis with Dissemination": 3, "Tuberculosis": 2},
        "Pulmonary nodules": {"Pulmonary Nocardiosis with Dissemination": 2, "Granulomatosis with Polyangiitis (GPA)": 2, "Metastatic Lung Cancer": 2, "Tuberculosis": 2},
        "Cutaneous lesions": {"Pulmonary Nocardiosis with Dissemination": 2, "Granulomatosis with Polyangiitis (GPA)": 1},
        "Neurological symptoms (confusion)": {"Pulmonary Nocardiosis with Dissemination": 2, "Metastatic Lung Cancer": 1},
        "Ineffective aminoglycoside therapy": {"Pulmonary Nocardiosis with Dissemination": 3}, # Nocardia is not treated with aminoglycosides
        "Productive cough & fever (acute infection)": {"Pulmonary Nocardiosis with Dissemination": 2, "Tuberculosis": 1},
        "History of ship building / smoking": {"Metastatic Lung Cancer": 2}
    }

    print("Synthesizing clinical findings to reach a diagnosis:")
    print("-" * 50)
    
    # Calculate scores for each diagnosis
    for finding, points in findings.items():
        for diagnosis, score in points.items():
            if diagnosis in diagnoses:
                diagnoses[diagnosis] += score

    # Find the diagnosis with the highest score
    most_likely_diagnosis = max(diagnoses, key=diagnoses.get)
    
    # Print the "equation" leading to the final diagnosis
    print(f"Final Diagnosis based on cumulative score: {most_likely_diagnosis}\n")
    print("Key Contributing Factors:")
    for finding, points in findings.items():
        if most_likely_diagnosis in points:
             # This loop prints out each number in the "equation"
             print(f"+ {points[most_likely_diagnosis]} for: {finding}")
             
    print("-" * 50)
    
    # Final explanation
    print(f"\nThe patient's presentation strongly suggests an opportunistic infection in an immunocompromised host. The combination of pulmonary nodules, neurological symptoms (confusion), and new cutaneous lesions after starting steroids is classic for a disseminated infection. Nocardia is a bacterium found in soil that fits this picture perfectly. Crucially, Nocardiosis does not respond to Aminoglycoside therapy, which explains the treatment failure and subsequent death from septic shock. While other conditions like GPA or cancer could be underlying, the proximate cause of death was most likely disseminated Nocardiosis.")


diagnose_patient()
<<<Pulmonary Nocardiosis with Dissemination>>>