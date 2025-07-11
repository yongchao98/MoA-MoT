def diagnose_patient():
    """
    This function evaluates potential diagnoses based on a patient's clinical and radiological findings
    using a weighted scoring system.
    """
    # Define the diagnoses and their corresponding letters
    diagnoses = {
        'A': 'Crohn\'s Disease',
        'B': 'Yersinia Colitis',
        'C': 'Ileocecal Tuberculosis',
        'D': 'Salmonella Enteritis',
        'E': 'C. difficile Colitis',
        'F': 'Cecal Volvulus',
        'G': 'Ischemic Colitis',
        'H': 'Pneumoperitoneum',
        'I': 'Ovarian Torsion',
        'J': 'Celiac Disease',
        'K': 'Gastrointestinal Lymphoma'
    }

    # Define key clinical findings and assign a weight based on diagnostic importance
    # The history of uveitis/arthritis is highly specific and thus weighted heavily.
    findings = {
        'ileocecal_inflammation': 4,
        'chronic_constitutional_symptoms': 3,
        'history_of_uveitis_and_arthritis': 5,
        'positive_inflammatory_markers': 2, # (Leukocytosis, +FOBT)
        'acute_on_chronic_pain': 2
    }

    # Define which findings are associated with each diagnosis
    diagnosis_features = {
        'Crohn\'s Disease': ['ileocecal_inflammation', 'chronic_constitutional_symptoms', 'history_of_uveitis_and_arthritis', 'positive_inflammatory_markers', 'acute_on_chronic_pain'],
        'Yersinia Colitis': ['ileocecal_inflammation', 'positive_inflammatory_markers'],
        'Ileocecal Tuberculosis': ['ileocecal_inflammation', 'chronic_constitutional_symptoms', 'positive_inflammatory_markers', 'acute_on_chronic_pain'],
        'Salmonella Enteritis': ['positive_inflammatory_markers'],
        'C. difficile Colitis': ['positive_inflammatory_markers'],
        'Cecal Volvulus': [],
        'Ischemic Colitis': ['positive_inflammatory_markers'],
        'Pneumoperitoneum': [],
        'Ovarian Torsion': [],
        'Celiac Disease': ['chronic_constitutional_symptoms'],
        'Gastrointestinal Lymphoma': ['ileocecal_inflammation', 'chronic_constitutional_symptoms', 'positive_inflammatory_markers']
    }

    # Calculate scores
    scores = {}
    for letter, name in diagnoses.items():
        score = 0
        if name in diagnosis_features:
            for feature in diagnosis_features[name]:
                score += findings[feature]
        scores[f"{letter}. {name}"] = score

    # Sort diagnoses by score in descending order
    sorted_diagnoses = sorted(scores.items(), key=lambda item: item[1], reverse=True)

    print("Likelihood Score for each Diagnosis based on Patient's Presentation:\n")
    for diagnosis, score in sorted_diagnoses:
        print(f"{diagnosis}: Score = {score}")

    most_likely = sorted_diagnoses[0]
    print(f"\nThe most likely diagnosis is {most_likely[0]} with a score of {most_likely[1]}.")

# Run the diagnostic evaluation
diagnose_patient()