import pandas as pd

def diagnose_patient_case():
    """
    Analyzes patient symptoms from a clinical case to suggest the most likely diagnosis.
    This function models a diagnostic process by scoring potential diseases
    based on the presence of key signs, symptoms, and risk factors.
    """

    # Key clinical clues extracted from the patient's case history
    patient_clues = {
        'age_over_60': True,
        'smoking_history': True,
        'asbestos_exposure_risk': True, # From shipbuilding
        'systemic_inflammation': True,   # Fatigue, polyarthralgia
        'pulmonary_nodules': True,
        'neurological_symptoms': True,   # Dizziness, confusion
        'hematological_symptoms': True,  # Bruising
        'cutaneous_lesions': True,
        'immunosuppressed_state': True, # Due to steroid treatment
        'terminal_sepsis': True,         # Died from septic shock
        'renal_involvement_signs': True, # Bruising, confusion can be signs of uremia
        'poor_response_to_antibiotics': True,
    }

    # Define diseases and their characteristic features for scoring
    disease_features = {
        'Granulomatosis with Polyangiitis (GPA)': [
            'systemic_inflammation', 'pulmonary_nodules', 'neurological_symptoms',
            'hematological_symptoms', 'cutaneous_lesions', 'immunosuppressed_state',
            'terminal_sepsis', 'renal_involvement_signs', 'poor_response_to_antibiotics'
        ],
        'Metastatic Lung Cancer': [
            'age_over_60', 'smoking_history', 'pulmonary_nodules',
            'neurological_symptoms', # Could be brain mets
            'hematological_symptoms' # Could be paraneoplastic
        ],
        'Miliary Tuberculosis': [
            'systemic_inflammation', 'pulmonary_nodules', 'terminal_sepsis', 'immunosuppressed_state'
        ],
        'Asbestosis/Mesothelioma': [
            'age_over_60', 'asbestos_exposure_risk', 'pulmonary_nodules'
        ]
    }

    # Calculate scores for each disease
    scores = {}
    explanation = {}
    for disease, features in disease_features.items():
        score = 0
        matched_features = []
        for feature in features:
            if patient_clues.get(feature, False):
                score += 1
                matched_features.append(feature)
        scores[disease] = score
        explanation[disease] = matched_features

    # Find the disease with the highest score
    most_likely_diagnosis = max(scores, key=scores.get)

    # Print the diagnostic reasoning process
    print("Clinical Case Analysis:")
    print("-----------------------")
    print("Based on the patient's presentation, a differential diagnosis was considered.")
    print("The following scores were calculated based on matching clinical features:\n")

    # Use pandas for a clean, tabular output
    results_df = pd.DataFrame(list(scores.items()), columns=['Diagnosis', 'Score']).sort_values(by='Score', ascending=False)
    print(results_df.to_string(index=False))

    print("\n-----------------------")
    print(f"Conclusion: The patient's constellation of symptoms, including systemic inflammation (fatigue, polyarthralgia), \n"
          f"multi-organ involvement (pulmonary nodules, confusion, bruising, cutaneous lesions), and progression to septic shock \n"
          f"in an immunosuppressed state, most strongly supports a diagnosis of:")
    print(f"\n>>> {most_likely_diagnosis}")

if __name__ == '__main__':
    diagnose_patient_case()
