def diagnose_patient():
    """
    This function simulates a diagnostic process by scoring potential diagnoses
    based on a patient's clinical, laboratory, and imaging findings.
    """

    # Patient's key findings with assigned weights for significance
    patient_findings = {
        'chronic_course': 2,        # Month-long duration
        'constitutional_symptoms': 2, # Fever, weight loss
        'ileocecal_inflammation': 3,  # CT finding of thickened ileum/cecum
        'uveitis_arthritis': 5,     # Highly specific extra-intestinal manifestations
        'rlq_pain': 1,              # Right lower quadrant pain
        'leukocytosis': 1,          # High WBC count
        'positive_fobt': 1,         # Blood in stool
        'age_60s': 1,               # Patient is 67
        'animal_exposure': 1        # Potential risk factor for some infections
    }

    # Characteristics of each diagnosis
    diagnoses = {
        'A. Crohn\'s Disease': ['chronic_course', 'constitutional_symptoms', 'ileocecal_inflammation', 'uveitis_arthritis', 'rlq_pain', 'leukocytosis', 'positive_fobt', 'age_60s'],
        'B. Yersinia Colitis': ['ileocecal_inflammation', 'rlq_pain', 'leukocytosis', 'animal_exposure'], # Typically acute, not chronic
        'C. Ileocecal Tuberculosis': ['chronic_course', 'constitutional_symptoms', 'ileocecal_inflammation', 'rlq_pain', 'positive_fobt'], # Great mimic, but lacks EIMs
        'D. Salmonella Enteritis': ['rlq_pain', 'leukocytosis'], # Acute
        'E. C. difficile Colitis': ['leukocytosis'], # Wrong location/chronicity
        'F. Cecal Volvulus': [], # Acute presentation, different CT
        'G. Ischemic Colitis': ['age_60s'], # Wrong location/chronicity
        'H. Pneumoperitoneum': [], # Not a primary diagnosis, not present
        'I. Ovarian Torsion': [], # Acute presentation, different pathology
        'J. Celiac Disease': ['chronic_course', 'constitutional_symptoms'], # Wrong location
        'K. Gastrointestinal Lymphoma': ['chronic_course', 'constitutional_symptoms', 'ileocecal_inflammation', 'positive_fobt', 'age_60s'] # Lacks EIMs
    }

    # Calculate scores
    scores = {}
    print("Calculating diagnostic scores based on patient findings:\n")
    for diagnosis, features in diagnoses.items():
        score = 0
        equation_parts = []
        for feature in features:
            if feature in patient_findings:
                points = patient_findings[feature]
                score += points
                equation_parts.append(str(points))
        
        scores[diagnosis] = score
        if not equation_parts:
            equation_parts.append("0")
            
        # Add a penalty for diagnoses that contradict the chronic course
        if "chronic_course" not in features and diagnosis in ['B. Yersinia Colitis', 'D. Salmonella Enteritis', 'F. Cecal Volvulus']:
             score -= 5 # Heavy penalty for being acute when presentation is chronic
             equation_parts.append("-5")

        scores[diagnosis] = score
        
        # Print the scoring for each diagnosis
        equation_str = " + ".join(equation_parts)
        print(f"- {diagnosis}: Score = {score} (Calculation: {equation_str})")


    # Find the most likely diagnosis
    most_likely_diagnosis = max(scores, key=scores.get)
    
    print("\n-------------------------------------------------")
    print(f"The most likely diagnosis is: {most_likely_diagnosis}")
    print("-------------------------------------------------\n")

# Run the diagnostic simulation
diagnose_patient()