def diagnose_patient():
    """
    Analyzes a clinical case to determine the most likely diagnosis
    by scoring potential diseases against key findings.
    """

    # Key findings from the clinical case description
    findings = {
        "Immunocompromised (Steroids)": 1,
        "Pulmonary Nodules": 1,
        "CNS Symptoms (Confusion)": 1,
        "Cutaneous Lesions": 1,
        "Aminoglycoside Therapy Failure": 1,
        "History (Smoking/Asbestos)": 1
    }

    # Scoring matrix for differential diagnoses.
    # A high positive score means the finding strongly supports the diagnosis.
    # A high negative score means the finding argues against the diagnosis.
    # 0 means the finding is neutral or irrelevant.
    disease_scores = {
        "Nocardiosis": {
            "Immunocompromised (Steroids)": 3,  # Classic risk factor
            "Pulmonary Nodules": 2,  # Classic presentation
            "CNS Symptoms (Confusion)": 2,  # Common site of dissemination
            "Cutaneous Lesions": 2,  # Common site of dissemination
            "Aminoglycoside Therapy Failure": 3,  # Characteristic resistance
            "History (Smoking/Asbestos)": 1  # Pre-existing lung damage is a risk factor
        },
        "Tuberculosis": {
            "Immunocompromised (Steroids)": 2,  # Major risk factor
            "Pulmonary Nodules": 2,  # Can present this way
            "CNS Symptoms (Confusion)": 1,  # Possible dissemination
            "Cutaneous Lesions": 1,  # Possible, but less common
            "Aminoglycoside Therapy Failure": -3, # Aminoglycosides are often used TO TREAT TB
            "History (Smoking/Asbestos)": 1 # Risk factor
        },
        "Lung Cancer": {
            "Immunocompromised (Steroids)": 0,  # Cancer itself is a risk, but this doesn't explain the acute infection
            "Pulmonary Nodules": 3,  # Hallmark sign
            "CNS Symptoms (Confusion)": 2,  # Common site for metastasis
            "Cutaneous Lesions": 1,  # Possible metastasis
            "Aminoglycoside Therapy Failure": 0, # Irrelevant for cancer itself
            "History (Smoking/Asbestos)": 3  # Major risk factors
        },
        "Wegener's Granulomatosis (GPA)": {
            "Immunocompromised (Steroids)": 0, # The treatment causes this, it's not part of the primary disease
            "Pulmonary Nodules": 2, # Classic feature
            "CNS Symptoms (Confusion)": 1, # Can occur
            "Cutaneous Lesions": 2, # Can occur
            "Aminoglycoside Therapy Failure": 0, # Irrelevant for an autoimmune disease
            "History (Smoking/Asbestos)": 0 # No direct link
        }
    }

    # Calculate final scores
    final_scores = {}
    for disease, scores in disease_scores.items():
        total_score = sum(scores[finding] for finding in findings)
        final_scores[disease] = total_score

    # Determine the most likely diagnosis
    most_likely_disease = max(final_scores, key=final_scores.get)

    print("Diagnostic Analysis based on Clinical Findings:\n")
    for disease, score in sorted(final_scores.items(), key=lambda item: item[1], reverse=True):
        print(f"- {disease}: Score = {score}")

    print("\n---\n")
    print(f"Conclusion: The constellation of symptoms, particularly the involvement of lungs, skin, and brain in an immunocompromised patient, combined with the failure of aminoglycoside therapy, most strongly suggests a specific diagnosis.")
    print(f"The highest scoring diagnosis is: {most_likely_disease}")

diagnose_patient()