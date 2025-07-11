def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring each option against the patient's findings.
    """
    # Step 1: Define the patient's clinical presentation
    patient = {
        'age': 1,
        'symptoms': {'hypertrophic scarring', 'erythema', 'spasticity'},
        'labs': {'anti-Mi-2': 'negative'}
    }

    # Step 2: Define features of potential diagnoses.
    # Scores: 2=Hallmark, 1=Common, 0=Not Associated, -1=Contradictory
    diagnoses_features = {
        'A. Ectropion': {
            'explanation': 'An outwardly turned eyelid; a physical sign, not a systemic disease to explain the full picture.',
            'hypertrophic scarring': 0, 'erythema': 0, 'spasticity': -1
        },
        'B. McArdle disease': {
            'explanation': 'A muscle metabolism disorder, typically presenting in adolescence with exercise intolerance, not with these skin findings or spasticity in infancy.',
            'hypertrophic scarring': -1, 'erythema': -1, 'spasticity': 0
        },
        'C. Dermatomyositis': {
            'explanation': 'Juvenile Dermatomyositis fits well. It causes inflammatory skin changes (erythema) which can ulcerate and lead to scarring, and muscle inflammation (myositis) which can manifest with weakness or spasticity. A negative anti-Mi-2 is common in the juvenile form.',
            'hypertrophic scarring': 1, 'erythema': 2, 'spasticity': 1
        },
        'D. McCune Albright syndrome': {
            'explanation': 'Characterized by fibrous dysplasia, cafe-au-lait spots, and endocrinopathy, which do not match the patient\'s presentation.',
            'hypertrophic scarring': -1, 'erythema': 0, 'spasticity': -1
        },
        'E. Cataracts': {
            'explanation': 'Clouding of the eye\'s lens; a finding, not a systemic diagnosis explaining the skin and muscle signs.',
            'hypertrophic scarring': -1, 'erythema': -1, 'spasticity': -1
        }
    }

    # Step 3 & 4: Calculate a score for each diagnosis and print the evaluation
    best_diagnosis = ''
    max_score = -float('inf')

    print("Evaluating diagnoses based on patient's findings...\n")

    for dx, features in diagnoses_features.items():
        score = 0
        equation = []

        # Score based on symptoms
        for symptom in patient['symptoms']:
            symptom_score = features.get(symptom, 0)
            score += symptom_score
            equation.append(f"{symptom.replace('_', ' ').title()}[{symptom_score}]")
            
        # Adjust score for lab results in specific context
        if dx == 'C. Dermatomyositis' and patient['labs']['anti-Mi-2'] == 'negative':
            # Negative anti-Mi-2 is expected in Juvenile Dermatomyositis, so it supports the diagnosis.
            lab_score = 1
            score += lab_score
            equation.append(f"Negative Anti-Mi-2 in Juvenile Patient[{lab_score}]")

        print(f"Diagnosis: {dx}")
        print(f"  - Rationale: {features['explanation']}")
        # Final formatting of the "equation" output
        final_equation_str = " + ".join(equation)
        print(f"  - Compatibility Score Equation: {final_equation_str} = {score}\n")

        if score > max_score:
            max_score = score
            best_diagnosis = dx

    # Step 5: Announce the most likely diagnosis
    print("------------------------------------------------------------")
    print(f"Conclusion: The most likely diagnosis is '{best_diagnosis}' with a score of {max_score}.")
    print("------------------------------------------------------------")

# Run the diagnostic analysis
solve_medical_case()