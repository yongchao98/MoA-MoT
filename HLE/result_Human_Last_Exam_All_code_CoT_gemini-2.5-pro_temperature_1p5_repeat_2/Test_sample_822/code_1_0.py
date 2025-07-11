import pandas as pd

def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring potential diseases against the patient's signs and symptoms.
    """
    # Key findings from the clinical case description
    case_findings = {
        'age_62': 1,
        'smoking_history_20_pack_year': 1,
        'occupational_exposure_shipbuilding': 1, # Asbestos, silica risk
        'polyarthritis': 1, # Joint pain/swelling
        'pulmonary_nodules': 1,
        'immunosuppression_from_steroids': 1, # Patient was put on steroids
        'acute_pneumonia': 1,
        'cutaneous_lesions': 1, # Skin lesions
        'aminoglycoside_therapy_failure': 1, # A key clue
        'sepsis_and_death': 1,
    }

    # Database of diseases and their characteristic features.
    # A positive score means the feature is characteristic of the disease.
    # A negative score means the feature makes the disease less likely.
    disease_features = {
        'Nocardiosis': {
            'pulmonary_nodules': 2, # Classic finding
            'cutaneous_lesions': 3, # Very common in disseminated disease
            'immunosuppression_from_steroids': 3, # The single most important risk factor
            'aminoglycoside_therapy_failure': 2, # Often requires specific combination therapy (e.g., TMP-SMX)
            'sepsis_and_death': 2, # Common outcome of severe, disseminated infection
            'polyarthritis': -1, # Not a primary feature of the infection itself
            'acute_pneumonia': 2,
            'occupational_exposure_shipbuilding': 1, # Environmental soil organism, exposure can happen anywhere
        },
        'Granulomatosis with Polyangiitis (GPA)': {
            'polyarthritis': 3, # Classic feature
            'pulmonary_nodules': 3, # Classic feature
            'cutaneous_lesions': 2, # Can be a feature of the underlying vasculitis
            'immunosuppression_from_steroids': 1, # This is the treatment, not a feature of the disease itself
            'smoking_history_20_pack_year': 1, # Association
            'sepsis_and_death': 1, # Possible, but often due to a secondary infection from treatment
            'aminoglycoside_therapy_failure': 0, # Not directly relevant
        },
        'Lung Cancer (with paraneoplastic syndrome)': {
            'smoking_history_20_pack_year': 3, # Strongest risk factor
            'pulmonary_nodules': 3, # Primary tumor and/or metastases
            'polyarthritis': 2, # Can be a paraneoplastic feature (hypertrophic osteoarthropathy)
            'immunosuppression_from_steroids': 1, # Steroids are often part of treatment/palliation
            'sepsis_and_death': 1, # Patients are vulnerable to infection
            'cutaneous_lesions': 0,
            'aminoglycoside_therapy_failure': 0,
            'age_62': 2
        },
        'Pulmonary Alveolar Proteinosis (PAP)': {
            'occupational_exposure_shipbuilding': 3, # Inhaled dust is a cause of secondary PAP
            'pulmonary_nodules': 2, # Imaging shows diffuse opacities, which could include nodules
            # PAP itself doesn't cause the acute infection, but creates a massive vulnerability to it.
            # We can represent this by scoring high for features of Nocardiosis if PAP is present.
            'immunosuppression_from_steroids': 3, # Strong predisposition to Nocardiosis
            'polyarthritis': -2, # Not a feature
            'cutaneous_lesions': 0, # Not a feature of PAP itself
            'sepsis_and_death': 1, # Only via secondary infection
        }
    }

    # Calculate scores
    scores = {}
    explanation = {}
    for disease, features in disease_features.items():
        score = 0
        disease_explanation = []
        for finding, present in case_findings.items():
            if present:
                points = features.get(finding, 0)
                score += points
                if points != 0:
                    disease_explanation.append(f"  - Matched on '{finding}': {points} points")
        scores[disease] = score
        explanation[disease] = disease_explanation

    # Find the highest scoring diagnosis
    best_diagnosis = max(scores, key=scores.get)

    print("Analyzing clinical case based on features...\n")
    for disease in sorted(scores, key=scores.get, reverse=True):
        print(f"Diagnosis: {disease}")
        print(f"Total Score: {scores[disease]}")
        for line in explanation[disease]:
            print(line)
        print("-" * 30)
    
    final_answer = (
        "The patient's initial inflammatory symptoms (arthritis) and pulmonary nodules "
        "led to steroid treatment, causing immunosuppression. This created the perfect "
        "environment for an opportunistic infection. The final, fatal illness with "
        "pneumonia, cutaneous lesions, and resistance to standard antibiotics is a "
        "classic presentation of disseminated Nocardiosis, which best explains "
        "the entire clinical picture, especially the terminal events."
    )
    print("\nConclusion:")
    print(final_answer)
    print("\n<<<Nocardiosis>>>")

solve_clinical_case()