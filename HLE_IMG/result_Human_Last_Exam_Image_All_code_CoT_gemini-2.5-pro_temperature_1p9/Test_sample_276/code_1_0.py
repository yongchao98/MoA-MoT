def diagnose_patient():
    """
    This function simulates a diagnostic scoring process based on patient data.
    """
    # Patient Data Points and their assigned weights for the scoring model
    clinical_features = {
        'age_gt_60': 1,
        'chronic_constitutional_symptoms': 3,
        'history_of_uveitis_and_arthritis': 3,
        'right_lower_quadrant_localization': 2,
        'leukocytosis_and_positive_fobt': 1,
        'ct_ileocecal_thickening': 3,
        'acute_presentation_negative_modifier': -3, # Penalty for chronic diseases presenting acutely
        'mismatch_ct_negative_modifier': -3, # Penalty for diagnoses that don't match CT
        'mismatch_history_negative_modifier': -2, # Penalty for diagnoses that don't match the history
    }

    # Scoring matrix for each diagnosis
    diagnoses_scores = {
        'A. Crohn\'s Disease':
            clinical_features['age_gt_60'] + \
            clinical_features['chronic_constitutional_symptoms'] + \
            clinical_features['history_of_uveitis_and_arthritis'] + \
            clinical_features['right_lower_quadrant_localization'] + \
            clinical_features['leukocytosis_and_positive_fobt'] + \
            clinical_features['ct_ileocecal_thickening'],

        'B. Yersinia Colitis':
            clinical_features['right_lower_quadrant_localization'] + \
            clinical_features['leukocytosis_and_positive_fobt'] + \
            clinical_features['mismatch_history_negative_modifier'], # History is chronic, not acute

        'C. Ileocecal Tuberculosis':
            clinical_features['age_gt_60'] + \
            clinical_features['chronic_constitutional_symptoms'] + \
            clinical_features['right_lower_quadrant_localization'] + \
            clinical_features['leukocytosis_and_positive_fobt'] + \
            clinical_features['ct_ileocecal_thickening'], # Lower score for uveitis/arthritis combo

        'K. Gastrointestinal Lymphoma':
            clinical_features['age_gt_60'] + \
            clinical_features['chronic_constitutional_symptoms'] + \
            clinical_features['right_lower_quadrant_localization'] + \
            clinical_features['leukocytosis_and_positive_fobt'] + \
            clinical_features['ct_ileocecal_thickening'] -1, # Uveitis/Arthritis less typical

        # Other diagnoses with significant negative scores
        'F. Cecal Volvulus': clinical_features['acute_presentation_negative_modifier'],
        'J. Celiac Disease': clinical_features['mismatch_ct_negative_modifier']
    }
    
    # Numbers for the equation as requested
    p1 = clinical_features['age_gt_60']
    p2 = clinical_features['chronic_constitutional_symptoms']
    p3 = clinical_features['history_of_uveitis_and_arthritis']
    p4 = clinical_features['right_lower_quadrant_localization']
    p5 = clinical_features['leukocytosis_and_positive_fobt']
    p6 = clinical_features['ct_ileocecal_thickening']

    print("Diagnostic Score Calculation:\n")
    for diagnosis, score in diagnoses_scores.items():
        print(f"{diagnosis}: Score = {score}")

    # Determine the most likely diagnosis
    most_likely_diagnosis = max(diagnoses_scores, key=diagnoses_scores.get)
    
    print("\n--- Conclusion ---")
    print(f"The most likely diagnosis based on the scoring is: {most_likely_diagnosis}")

    # Final equation for the top diagnosis
    final_score = diagnoses_scores[most_likely_diagnosis]
    print(f"\nFinal scoring equation for Crohn's Disease:")
    print(f"Score = {p1} (Age) + {p2} (Symptoms) + {p3} (History) + {p4} (Exam) + {p5} (Labs) + {p6} (CT) = {final_score}")


diagnose_patient()