def solve_diagnosis():
    """
    Analyzes a patient's clinical presentation to determine the most likely diagnosis
    from a given set of options by scoring each option against the patient's data.
    """
    patient_profile = {
        'age': 1,
        'symptoms': ['hypertrophic scarring', 'erythema', 'spasticity (or severe contractures)'],
        'labs': 'anti-Mi-2 negative'
    }

    diagnoses = {
        'A. Ectropion': {
            'features': ['eyelid turning outward'],
            'mismatches': ['hypertrophic scarring', 'spasticity (or severe contractures)'],
            'score_details': "Fails to explain systemic symptoms."
        },
        'B. McArdle disease': {
            'features': ['exercise intolerance', 'muscle cramps'],
            'mismatches': ['hypertropehic scarring', 'erythema', 'infant presentation'],
            'score_details': "Wrong age of onset and does not fit skin findings."
        },
        'C. Dermatomyositis': {
            'features': ['erythema', 'muscle involvement (can lead to contractures/spasticity)',
                         'scarring (from skin ulceration)', 'anti-Mi-2 negative (common in JDM)'],
            'mismatches': [],
            'score_details': "Consistent with Juvenile Dermatomyositis (JDM)."
        },
        'D. McCune Albright syndrome': {
            'features': ['cafe-au-lait spots', 'fibrous dysplasia', 'endocrinopathy'],
            'mismatches': ['erythema', 'spasticity (or severe contractures)', 'hypertrophic scarring'],
            'score_details': "Presents with a different set of classic signs."
        },
        'E. Cataracts': {
            'features': ['clouding of eye lens'],
            'mismatches': ['hypertrophic scarring', 'erythema', 'spasticity (or severe contractures)'],
            'score_details': "A single symptom, not a unifying diagnosis."
        }
    }

    best_match = None
    highest_score = -1
    final_equation_str = ""

    print("Analyzing patient profile against possible diagnoses:\n")

    for diagnosis, data in diagnoses.items():
        score = 0
        equation_parts = []
        
        # Score based on matching symptoms
        for symptom in patient_profile['symptoms']:
            if any(symptom_keyword in feature for feature in data['features'] for symptom_keyword in symptom.split()):
                score += 1
                equation_parts.append(f"1 ({symptom.split('(')[0].strip()})")

        # Score based on matching lab results, if relevant
        if 'anti-Mi-2 negative' in data.get('features', []):
            score += 1
            equation_parts.append("1 (Anti-Mi-2 Negative)")

        # Penalize for mismatches
        for mismatch in data.get('mismatches', []):
            score -= 1
        
        print(f"Diagnosis: {diagnosis}")
        print(f"Rationale: {data['score_details']}")
        print(f"Compatibility Score: {score}\n")

        if score > highest_score:
            highest_score = score
            best_match = diagnosis
            final_equation_str = f"Final Score for {best_match}: {' + '.join(equation_parts)} = {score}"


    print("--- Conclusion ---")
    print(f"The most likely diagnosis is {best_match} with a score of {highest_score}.")
    print("The final calculation for the best match is:")
    print(final_equation_str)

solve_diagnosis()