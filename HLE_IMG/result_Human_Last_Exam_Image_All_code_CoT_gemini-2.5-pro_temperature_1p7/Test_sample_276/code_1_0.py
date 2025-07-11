def diagnose_patient():
    """
    This function models the diagnostic process by scoring different potential diagnoses
    based on the patient's clinical and radiological findings.
    """

    # Define the patient's key findings from the case description
    findings = {
        'chronic_symptoms': {'value': True, 'score_weight': 2, 'desc': "chronic symptoms (1 month)"},
        'constitutional_symptoms': {'value': True, 'score_weight': 2, 'desc': "weight loss/fever"},
        'location_rlq': {'value': True, 'score_weight': 2, 'desc': "RLQ pain/tenderness"},
        'leukocytosis': {'value': True, 'score_weight': 1, 'desc': "WBC 13,000"},
        'fobt_positive': {'value': True, 'score_weight': 1, 'desc': "positive FOBT"},
        'extra_intestinal': {'value': True, 'score_weight': 3, 'desc': "uveitis/arthritis"},
        'ct_ileocecal_thickening': {'value': True, 'score_weight': 2, 'desc': "ileocecal thickening on CT"},
        'acute_presentation': {'value': False, 'score_weight': -2, 'desc': "strictly acute onset"}
    }

    # Define how each diagnosis aligns with the findings.
    # Score: + for consistent, - for inconsistent, 0 for neutral.
    diagnosis_matrix = {
        "Crohn's Disease": {
            'chronic_symptoms': 1, 'constitutional_symptoms': 1, 'location_rlq': 1,
            'leukocytosis': 1, 'fobt_positive': 1, 'extra_intestinal': 1,
            'ct_ileocecal_thickening': 1, 'acute_presentation': 0
        },
        "Yersinia Colitis": {
            'chronic_symptoms': -1, 'constitutional_symptoms': 0, 'location_rlq': 1,
            'leukocytosis': 1, 'fobt_positive': 0, 'extra_intestinal': -0.5,
            'ct_ileocecal_thickening': 1, 'acute_presentation': 1
        },
        "Ileocecal Tuberculosis": {
            'chronic_symptoms': 1, 'constitutional_symptoms': 1, 'location_rlq': 1,
            'leukocytosis': 1, 'fobt_positive': 0, 'extra_intestinal': -0.5,
            'ct_ileocecal_thickening': 1, 'acute_presentation': 0
        },
        "Cecal Volvulus": {
            'chronic_symptoms': -1, 'constitutional_symptoms': -1, 'location_rlq': 1,
            'leukocytosis': 0, 'fobt_positive': 0, 'extra_intestinal': -1,
            'ct_ileocecal_thickening': -1, 'acute_presentation': 1
        },
        "GI Lymphoma": {
            'chronic_symptoms': 1, 'constitutional_symptoms': 1, 'location_rlq': 1,
            'leukocytosis': 0, 'fobt_positive': 1, 'extra_intestinal': -1,
            'ct_ileocecal_thickening': 1, 'acute_presentation': 0
        }
    }

    # Calculate scores
    results = {}
    for diagnosis, multipliers in diagnosis_matrix.items():
        total_score = 0
        equation_parts = []
        for finding_name, finding_data in findings.items():
            if finding_data['value']: # Only consider positive findings for this model
                multiplier = multipliers.get(finding_name, 0)
                score = finding_data['score_weight'] * multiplier
                if score != 0:
                    total_score += score
                    equation_parts.append(f"{score} ({finding_data['desc']})")
        results[diagnosis] = {'score': total_score, 'equation': " + ".join(equation_parts)}

    # Find the best diagnosis
    best_diagnosis = max(results, key=lambda k: results[k]['score'])
    
    print("Diagnostic Score Calculation for the Most Likely Diagnosis:\n")
    print(f"Diagnosis: {best_diagnosis}")
    equation = results[best_diagnosis]['equation'].replace("+ -", "- ")
    final_score = results[best_diagnosis]['score']
    print(f"Score Calculation = {equation} = {final_score}")

    print("\nConclusion: The combination of chronic symptoms, extraintestinal manifestations (uveitis, arthritis), and classic ileocecal inflammation on CT makes Crohn's Disease the most likely diagnosis.")


diagnose_patient()