import pandas as pd

def diagnose_skin_condition():
    """
    Analyzes patient data to provide a differential diagnosis for skin lesions
    based on a scoring model.
    """
    # 1. Key clinical features extracted from the patient case
    patient_profile = {
        'obesity_bmi_39': 39,
        'smoker': True,
        'location_intertriginous': True,  # Axillary, inframammary, inguinal
        'lesion_purulent_nodules': True,
        'lesion_erythematous_plaques': True,
        'lesion_large_bullae_or_abscesses': True # "Large bullae" often describes large, fluctuant abscesses in HS
    }

    # 2. Scoring matrix: How strongly each feature points to a diagnosis.
    # High score = strongly suggests; negative score = argues against.
    scoring_rules = {
        'Hidradenitis Suppurativa': {
            'lesion_purulent_nodules': 30,
            'location_intertriginous': 20,
            'smoker': 15,
            'obesity_bmi_39': 15,
            'lesion_large_bullae_or_abscesses': 10,
            'lesion_erythematous_plaques': 5, # Represents surrounding inflammation
        },
        'Psoriasis (Inverse)': {
            'lesion_erythematous_plaques': 25,
            'location_intertriginous': 20,
            'obesity_bmi_39': 10,
            'smoker': 5,
            'lesion_purulent_nodules': -20, # Not typical for inverse psoriasis
            'lesion_large_bullae_or_abscesses': -15, # Not typical
        },
        'Allergic Contact Dermatitis': {
            'lesion_erythematous_plaques': 10,
            'location_intertriginous': 5,
            'lesion_purulent_nodules': -30, # Strong contraindication
            'lesion_large_bullae_or_abscesses': -20, # Not characteristic unless severe/infected
            'smoker': 0,
            'obesity_bmi_39': 0,
        },
        'Malignant Intertrigo': {
            'location_intertriginous': 10,
            'lesion_erythematous_plaques': 10,
            'lesion_purulent_nodules': -15, # Unlikely presentation
            'lesion_large_bullae_or_abscesses': -10, # Unlikely presentation
            'smoker': 0,
            'obesity_bmi_39': 0,
        },
         'Atopic Dermatitis': {
            'location_intertriginous': 5,
            'lesion_erythematous_plaques': 10,
            'lesion_purulent_nodules': -25, # Not typical
            'lesion_large_bullae_or_abscesses': -20, # Not typical
            'smoker': 0,
            'obesity_bmi_39': 0,
        }
    }

    # 3. Calculate and display scores
    results = {}
    print("Calculating Likelihood Scores...\n")
    
    diagnoses_map = {
        'Hidradenitis Suppurativa': 'C. Hidradenitis Suppurativa',
        'Psoriasis (Inverse)': 'E. Psoriasis',
        'Allergic Contact Dermatitis': 'B. Allergic contact dermatitis',
        'Malignant Intertrigo': 'A. Malignant Intertrigo',
        'Atopic Dermatitis': 'D. Atopic dermatitis'
    }

    for diagnosis, rules in scoring_rules.items():
        total_score = 0
        equation_str = f"Score({diagnoses_map[diagnosis]}) = "
        terms = []
        for feature, present in patient_profile.items():
            if present:
                score = rules.get(feature, 0)
                if score != 0:
                    terms.append(f"{score} ({feature})")
                total_score += score
        
        equation_str += " + ".join(terms).replace("+ -", "- ")
        print(f"Equation: {equation_str}")
        print(f"Final Score: {total_score}\n")
        results[diagnoses_map[diagnosis]] = total_score

    # 4. Determine the final diagnosis
    most_likely_diagnosis = max(results, key=results.get)

    print("-" * 30)
    print(f"Most Likely Diagnosis: {most_likely_diagnosis}")
    print(f"Confidence Score: {results[most_likely_diagnosis]}")
    print("-" * 30)


if __name__ == '__main__':
    diagnose_skin_condition()