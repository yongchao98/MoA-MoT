def diagnose_patient():
    """
    This script analyzes the clinical findings from the case study to determine the most likely diagnosis.
    It uses a scoring system to weigh the importance of each finding.
    """
    # Define potential diagnoses and their initial scores.
    scores = {
        "Lung Cancer": 0,
        "Granulomatosis with Polyangiitis (GPA)": 0,
        "Tuberculosis": 0,
        "Nocardiosis": 0
    }

    # Define clinical findings and their associated point values for each diagnosis.
    # Highly specific findings are given more weight.
    findings = {
        "Asbestos/Smoking History (Risk for malignancy)": {"Lung Cancer": 2},
        "Multiple Pulmonary Nodules": {"Lung Cancer": 2, "GPA": 2, "Nocardiosis": 2, "Tuberculosis": 2},
        "Immunocompromised (due to steroids/illness)": {"Nocardiosis": 3, "Tuberculosis": 2},
        "Cutaneous Lesions (Infectious)": {"Nocardiosis": 3, "GPA": 1},
        "CNS Symptoms (Confusion)": {"Nocardiosis": 3, "Lung Cancer": 1},
        "Aminoglycoside Therapy Ineffective": {"Nocardiosis": 5} # This is a highly specific clue
    }

    print("Diagnostic Analysis based on Clinical Findings:")
    print("="*45)

    # Calculate the scores based on the findings. This represents the "equation".
    for finding, disease_weights in findings.items():
        print(f"\nFinding: '{finding}'")
        for disease, weight in disease_weights.items():
            # This represents one term in the final equation for the score
            print(f"  - Adding {weight} to {disease}'s score.")
            scores[disease] += weight

    print("\n" + "="*45)
    print("Final Score Calculation:")
    # The final "equation" is the sum of all points for each disease.
    for disease, score in scores.items():
        # We print each part of the final equation.
        print(f"Total Score for {disease}: {score}")
    print("="*45)

    # Determine the most likely diagnosis by finding the highest score.
    most_likely_disease = max(scores, key=scores.get)

    print(f"\nConclusion:")
    print(f"The combination of an immunocompromised state with a triad of pulmonary, cutaneous,")
    print(f"and CNS symptoms, along with resistance to aminoglycoside therapy,")
    print(f"points most strongly to one diagnosis.")
    print(f"\nThe most likely disease is: {most_likely_disease}")

diagnose_patient()