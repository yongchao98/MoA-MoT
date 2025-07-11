import json

def evaluate_patient_case():
    """
    Analyzes a clinical case by scoring potential diagnoses against patient history.
    """
    patient_features = {
        'occupation': 'metal smelting (30 years)',
        'chronic_conditions': ['hypertension'],
        'initial_symptoms': ['mania-like symptoms'],
        'final_symptom': 'sexual dysfunction',
        'family_history': 'mood disorders',
        'likely_treatment': 'lithium'
    }

    diagnoses = {
        'A': {
            'name': 'Lithium induced hypothyroidism',
            'explains': ['final_symptom', 'likely_treatment', 'family_history']
        },
        'B': {
            'name': 'Arsenic induced Renal Dysfunction',
            'explains': ['occupation']
        },
        'C': {
            'name': 'Mercury induced Renal Dysfunction',
            'explains': ['occupation']
        },
        'D': {
            'name': 'Lead induced Sexual dysfunction',
            'explains': ['occupation', 'chronic_conditions', 'initial_symptoms', 'final_symptom']
        },
        'E': {
            'name': 'Manganese induced Renal Dysfunction',
            'explains': ['occupation']
        }
    }

    print("Analyzing the patient case by scoring potential root causes...")
    print("-" * 30)

    scores = {}
    for key, diagnosis in diagnoses.items():
        score = 0
        reasons = []
        for feature in diagnosis['explains']:
            if feature in patient_features:
                score += 1
                reasons.append(f"explains {feature} ('{patient_features[feature][0]}')")
            # Special case for 'likely_treatment' which is an inference
            elif feature == 'likely_treatment' or feature == 'family_history':
                score += 1
                reasons.append(f"is consistent with the inferred '{feature}'")
        
        # Unifying diagnosis bonus: a diagnosis that links occupation to symptoms gets an extra point
        if 'occupation' in diagnosis['explains'] and ('initial_symptoms' in diagnosis['explains'] or 'final_symptom' in diagnosis['explains']):
            score += 1
            reasons.append("provides a strong unifying link between occupation and symptoms")

        scores[key] = {'score': score, 'name': diagnosis['name'], 'reasoning': reasons}

    # Find the best diagnosis
    best_diagnosis_key = max(scores, key=lambda k: scores[k]['score'])
    best_diagnosis = scores[best_diagnosis_key]
    
    print("Scoring complete. The most likely root cause is:\n")
    print(f"Choice {best_diagnosis_key}: {best_diagnosis['name']}")
    print(f"Score: {best_diagnosis['score']}")
    print("Reasoning:")
    for reason in best_diagnosis['reasoning']:
        print(f"- {reason}")
        
    print("\nConclusion: Lead exposure provides the most comprehensive explanation, linking the patient's long-term occupation to his hypertension, initial neuropsychiatric symptoms, and eventual sexual dysfunction.")

evaluate_patient_case()