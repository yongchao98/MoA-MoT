def solve_clinical_case():
    """
    This script evaluates the clinical case by assigning scores to potential diagnoses
    based on the presence and significance of key findings.
    """

    print("Analyzing the clinical case based on key findings...\n")

    # These weights are assigned based on the clinical significance of each finding for the diagnosis.
    findings = {
        'Pulmonary Nodules': {'nocardiosis': 2, 'gpa': 3, 'lung_cancer': 4},
        'Immunosuppression (Steroids)': {'nocardiosis': 5, 'gpa': -2, 'lung_cancer': 0},
        'Cutaneous Lesions': {'nocardiosis': 4, 'gpa': 2, 'lung_cancer': 1},
        'Neurological Symptoms (Confusion)': {'nocardiosis': 4, 'gpa': 1, 'lung_cancer': 2},
        'Failed Aminoglycoside Therapy': {'nocardiosis': 5, 'gpa': 0, 'lung_cancer': 0},
        'Smoking & Asbestos History': {'nocardiosis': 0, 'gpa': 1, 'lung_cancer': 4},
        'Polyarthritis': {'nocardiosis': 0, 'gpa': 3, 'lung_cancer': 1},
        'Acute Sepsis Picture': {'nocardiosis': 3, 'gpa': 1, 'lung_cancer': 1}
    }

    # Note: GPA score is negative for steroids because GPA is the *reason* for the treatment,
    # not a consequence of it. We are diagnosing the final acute illness.

    diagnoses = ['nocardiosis', 'gpa', 'lung_cancer']
    scores = {dx: 0 for dx in diagnoses}
    equations = {dx: [] for dx in diagnoses}

    for finding, weights in findings.items():
        for dx in diagnoses:
            score_value = weights[dx]
            if score_value != 0:
                scores[dx] += score_value
                equations[dx].append(str(score_value))

    print("--- Diagnostic Score Calculation ---\n")

    # Nocardiosis Score
    nocardiosis_eq_str = " + ".join(equations['nocardiosis'])
    print(f"Disseminated Nocardiosis Score Equation: {nocardiosis_eq_str}")
    print(f"This score is based on the weights for: Pulmonary Nodules, Immunosuppression, Cutaneous Lesions, Neurological Symptoms, Failed Specific Therapy, and Acute Sepsis.")
    final_nocardiosis_score = scores['nocardiosis']
    print(f"Final Score for Disseminated Nocardiosis = {final_nocardiosis_score}\n")

    # GPA Score
    gpa_eq_str = " + ".join(equations['gpa']).replace("+ -", "- ")
    print(f"Granulomatosis with Polyangiitis (GPA) Score Equation: {gpa_eq_str}")
    print(f"This score is based on the weights for: Pulmonary Nodules, Steroid Use (as a pre-existing factor), Cutaneous Lesions, Neurological Symptoms, Exposure History, Polyarthritis, and Acute Sepsis.")
    final_gpa_score = scores['gpa']
    print(f"Final Score for GPA = {final_gpa_score}\n")

    # Lung Cancer Score
    lung_cancer_eq_str = " + ".join(equations['lung_cancer'])
    print(f"Lung Cancer Score Equation: {lung_cancer_eq_str}")
    print(f"This score is based on the weights for: Pulmonary Nodules, Neurological Symptoms, Exposure History, Polyarthritis (paraneoplastic), and Acute Sepsis.")
    final_lung_cancer_score = scores['lung_cancer']
    print(f"Final Score for Lung Cancer = {final_lung_cancer_score}\n")

    print("--- Conclusion ---")
    print(f"The highest score ({final_nocardiosis_score}) points towards Disseminated Nocardiosis as the direct cause of the patient's final, fatal illness. The infection was likely opportunistic due to steroid-induced immunosuppression, which was prescribed for an underlying condition (possibly GPA or malignancy).")


solve_clinical_case()
<<<Disseminated Nocardiosis>>>