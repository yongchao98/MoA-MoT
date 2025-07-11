def analyze_patient_data():
    """
    Analyzes patient data to suggest a possible diagnosis based on a simple scoring model.
    This is for illustrative purposes and is not a substitute for professional medical advice.
    """
    # --- Patient Findings ---
    # We represent each key finding and its assigned weight for the scoring model.
    # The vertebral mass is a highly specific and significant finding pointing towards metastasis.
    findings = {
        'vertebral_mass': {'present': True, 'weight': 10},
        'chronic_cough': {'present': True, 'weight': 3},
        'dyspnea': {'present': True, 'weight': 3},
        'copd_history': {'present': True, 'weight': 2},  # A significant risk factor
        'acid_reflux': {'present': True, 'weight': 1},
        'high_creatinine': {'present': True, 'weight': 1}
    }

    # --- Differential Diagnoses and their associated findings ---
    # We define which findings are typically associated with each condition.
    diagnoses = {
        'A. Aspiration pneumonitis': ['chronic_cough', 'dyspnea', 'acid_reflux'],
        'B. Aspiration pneumonia': ['chronic_cough', 'dyspnea', 'acid_reflux'],
        'C. Achalasia': ['acid_reflux'],
        'D. Adenocarcinoma': ['vertebral_mass', 'chronic_cough', 'dyspnea', 'copd_history', 'high_creatinine'],
        'E. COPD': ['chronic_cough', 'dyspnea']
    }

    scores = {}
    best_diagnosis = ""
    max_score = -1
    calculation_string = ""

    print("Analyzing patient findings against potential diagnoses...\n")

    # Calculate score for each diagnosis
    for diagnosis, associated_findings in diagnoses.items():
        current_score = 0
        current_calc_parts = []
        for finding_name in associated_findings:
            if findings[finding_name]['present']:
                weight = findings[finding_name]['weight']
                current_score += weight
                current_calc_parts.append(f"{weight} ({finding_name})")

        scores[diagnosis] = current_score

        if current_score > max_score:
            max_score = current_score
            best_diagnosis = diagnosis
            calculation_string = " + ".join(current_calc_parts) + f" = {max_score}"

    print(f"The most likely diagnosis based on the scoring model is: {best_diagnosis}")
    print("\nThis diagnosis provides the most comprehensive explanation for the patient's key findings, especially the vertebral mass, which is indicative of metastatic disease.")
    print("\nThe scoring was calculated as follows for the highest-scoring diagnosis:")
    print(f"Final Equation: {calculation_string}")

analyze_patient_data()