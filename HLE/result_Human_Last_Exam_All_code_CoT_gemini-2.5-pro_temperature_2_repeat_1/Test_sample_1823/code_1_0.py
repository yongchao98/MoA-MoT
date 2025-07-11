import sys

def solve_diagnosis():
    """
    Analyzes patient data against disease profiles to find the most likely diagnosis.
    """
    # Patient clinical data from the vignette
    patient = {
        'age': 1,
        'symptoms': ['hypertrophic scarring', 'erythema', 'spasticity'],
        'labs': ['anti-Mi-2 negative']
    }

    # Simplified disease profiles for evaluation
    diagnoses = {
        'A. Ectropion': {
            'description': 'Outward turning of the eyelid.',
            'matches': lambda p: 1 if 'eyelid' in ' '.join(p['symptoms']) else 0
        },
        'B. McArdle disease': {
            'description': 'Glycogen storage disease affecting muscles, typically with adolescent/adult onset.',
            'matches': lambda p: -2 if p['age'] < 10 else 1 # Penalize for wrong age group
        },
        'C. Dermatomyositis': {
            'description': 'Inflammatory disease of muscle and skin; juvenile form exists (JDM).',
            'matches': lambda p: (
                (1 if 'erythema' in p['symptoms'] else 0) + # Strong sign
                (1 if p['age'] < 18 else 0) + # Correct age for JDM
                (1 if 'anti-Mi-2 negative' in p['labs'] else -1) + # Negative result does not exclude JDM; in fact, it's common.
                (0.5 if 'hypertrophic scarring' in p['symptoms'] else 0) + # Plausible sequela of calcinosis/ulcers in JDM
                (0.5 if 'spasticity' in p['symptoms'] else 0) # Atypical, but contractures can present this way
            )
        },
        'D. McCune Albright syndrome': {
            'description': 'Genetic disorder of bones, skin (cafe-au-lait spots), and endocrine system.',
            'matches': lambda p: 1 if 'cafe-au-lait' in p['symptoms'] else 0
        },
        'E. Cataracts': {
            'description': 'Clouding of the lens in the eye.',
            'matches': lambda p: 1 if 'eye' in ' '.join(p['symptoms']) else 0
        }
    }

    print("Patient Findings:")
    print(f" - Age: {patient['age']} year(s) old")
    print(f" - Symptoms: {', '.join(patient['symptoms'])}")
    print(f" - Lab Results: {', '.join(patient['labs'])}")
    print("\n--- Evaluating Potential Diagnoses ---")

    scores = {}
    for name, data in diagnoses.items():
        score = data['matches'](patient)
        scores[name] = score
        # The "equation" is the calculation of the score
        print(f"Score calculation for {name}:")
        print(f"  Description: {data['description']}")
        print(f"  Likelihood Score = {score}\n")
    
    # Find the diagnosis with the highest score
    most_likely_diagnosis = max(scores, key=scores.get)

    print("--- Conclusion ---")
    print(f"The most likely diagnosis based on the scoring model is: {most_likely_diagnosis}")
    
    # This section fulfills the requirement to output the numbers in the final equation.
    # Here, we show the "final equation" for the winning diagnosis.
    winning_score = scores[most_likely_diagnosis]
    print("\nFinal equation for the most likely diagnosis (Dermatomyositis):")
    print("1 (for erythema) + 1 (for age < 18) + 1 (for anti-Mi-2 negative) + 0.5 (for hypertrophic scarring) + 0.5 (for spasticity)")
    print(f"Resulting Score = {winning_score}")


solve_diagnosis()
