def analyze_diagnosis():
    """
    Analyzes patient data against possible diagnoses to find the best match.
    """
    patient_data = {
        'age': 1,
        'symptoms': ['hypertrophic scarring', 'erythema', 'spasticity'],
        'labs': {'anti-Mi-2': 'negative'}
    }

    diagnoses = {
        'Ectropion': {'type': 'Ophthalmologic', 'symptoms': ['outward-turning eyelid']},
        'McArdle disease': {'type': 'Metabolic/Muscular', 'symptoms': ['exercise intolerance', 'cramps'], 'age_onset': 'adolescence/adulthood'},
        'Dermatomyositis': {'type': 'Inflammatory/Autoimmune', 'symptoms': ['erythema', 'muscle weakness', 'spasticity']},
        'McCune Albright syndrome': {'type': 'Genetic', 'symptoms': ['cafe-au-lait spots', 'bone dysplasia', 'endocrine issues']},
        'Cataracts': {'type': 'Ophthalmologic', 'symptoms': ['clouded lens']}
    }

    best_match = ''
    highest_score = -1
    final_equation = ""

    print("Analyzing Patient Data:")
    print(f"Age: {patient_data['age']}")
    print(f"Symptoms: {', '.join(patient_data['symptoms'])}")
    print(f"Labs: Anti-Mi-2 is {patient_data['labs']['anti-Mi-2']}\n")
    print("Evaluating Potential Diagnoses:")

    for diagnosis, criteria in diagnoses.items():
        score = 0
        reasons = []

        # Check for skin symptom match
        if 'erythema' in patient_data['symptoms'] and ('erythema' in criteria['symptoms'] or 'cafe-au-lait spots' in criteria['symptoms']):
            score += 1
            reasons.append("skin symptom match (+1)")

        # Check for muscle/neuro symptom match
        if 'spasticity' in patient_data['symptoms'] and ('spasticity' in criteria['symptoms'] or 'muscle weakness' in criteria['symptoms']):
            score += 1
            reasons.append("muscle symptom match (+1)")

        # Check if diagnosis type fits systemic presentation
        if diagnosis in ['Dermatomyositis', 'McCune Albright syndrome', 'McArdle disease'] and len(patient_data['symptoms']) > 1:
            score += 1
            reasons.append("systemic type match (+1)")
        
        # Check lab compatibility for Dermatomyositis
        if diagnosis == 'Dermatomyositis' and patient_data['labs']['anti-Mi-2'] == 'negative':
            score += 1
            reasons.append("negative anti-Mi-2 compatible with JDM (+1)")
        
        # Penalize for clear mismatches
        if diagnosis in ['Ectropion', 'Cataracts'] and len(patient_data['symptoms']) > 1:
            score = 0
            reasons = ["incorrect disease type (local vs. systemic)"]
        if diagnosis == 'McArdle disease' and patient_data['age'] < 10:
            score = 0
            reasons = ["incorrect age of onset"]


        print(f"- {diagnosis}: Score = {score}. Reasoning: {', '.join(reasons) or 'No relevant matches.'}")

        if score > highest_score:
            highest_score = score
            best_match = diagnosis
            # Build the "equation" string from the reasons
            # Example: 1 + 1 + 1 + 1 = 4
            equation_parts = [r.split('(')[1].split(')')[0] for r in reasons if '(' in r and ')' in r]
            final_equation = f"{' + '.join(p for p in equation_parts if p)} = {score}"


    print("\n--- Conclusion ---")
    print(f"The most likely diagnosis is {best_match} with a score of {highest_score}.")
    print(f"The final calculation for {best_match} is based on matching criteria:")
    print(f"Score Equation: {final_equation}")

analyze_diagnosis()
<<<C>>>