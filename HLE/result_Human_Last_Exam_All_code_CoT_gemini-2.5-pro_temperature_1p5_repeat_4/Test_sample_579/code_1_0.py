import pandas as pd

def diagnose_skin_condition():
    """
    Analyzes a clinical vignette to determine the most likely dermatological diagnosis
    using a weighted scoring system.
    """
    # Key findings from the clinical case
    findings = {
        'location_intertriginous': {'weight': 3, 'description': "Classic location (axillary, inframammary, inguinal)"},
        'lesion_purulent_nodules': {'weight': 5, 'description': "Hallmark finding of purulent nodules"},
        'lesion_plaques': {'weight': 2, 'description': "Presence of erythematous plaques"},
        'lesion_bullae': {'weight': 1, 'description': "Presence of large bullae"},
        'risk_factors': {'weight': 3, 'description': "Associated risk factors (obesity, smoking)"}
    }

    # Scoring matrix for each diagnosis based on how well it matches the findings (0=No Match, 1=Weak, 2=Moderate, 3=Strong)
    diagnoses_scores = {
        'A. Malignant Intertrigo':    {'location_intertriginous': 2, 'lesion_purulent_nodules': 0, 'lesion_plaques': 2, 'lesion_bullae': 0, 'risk_factors': 0},
        'B. Allergic contact dermatitis': {'location_intertriginous': 1, 'lesion_purulent_nodules': 0, 'lesion_plaques': 2, 'lesion_bullae': 1, 'risk_factors': 0},
        'C. Hidradenitis Suppurativa': {'location_intertriginous': 3, 'lesion_purulent_nodules': 3, 'lesion_plaques': 1, 'lesion_bullae': 1, 'risk_factors': 3},
        'D. Atopic dermatitis':       {'location_intertriginous': 1, 'lesion_purulent_nodules': 0, 'lesion_plaques': 2, 'lesion_bullae': 0, 'risk_factors': 0},
        'E. Psoriasis':               {'location_intertriginous': 3, 'lesion_purulent_nodules': 0, 'lesion_plaques': 3, 'lesion_bullae': 0, 'risk_factors': 1}
    }

    # Calculate final scores and prepare for display
    results = []
    for diagnosis, scores in diagnoses_scores.items():
        total_score = 0
        equation_parts = []
        for finding, score_value in scores.items():
            weighted_score = score_value * findings[finding]['weight']
            total_score += weighted_score
            equation_parts.append(str(weighted_score))
        
        results.append({
            'Diagnosis': diagnosis,
            'Total Score': total_score,
            'Equation': ' + '.join(equation_parts)
        })

    # Create a DataFrame for nice printing
    df = pd.DataFrame(results)
    df = df.sort_values(by='Total Score', ascending=False).reset_index(drop=True)

    # Print the analysis
    print("--- Diagnostic Scoring Analysis ---")
    print(df)
    
    # Get the top diagnosis
    top_diagnosis = df.iloc[0]
    best_match_name = top_diagnosis['Diagnosis']
    
    # Deconstruct the winning score to show the 'equation' as requested
    winning_scores = diagnoses_scores[best_match_name]
    
    print("\n--- Conclusion ---")
    print(f"The most likely diagnosis is {best_match_name} based on the scoring.")
    print("\nThe final score is based on how well the diagnosis matches key findings:")
    
    equation_str_parts = []
    final_score = 0
    for finding_key, score_value in winning_scores.items():
        weight = findings[finding_key]['weight']
        description = findings[finding_key]['description']
        calculated_value = score_value * weight
        final_score += calculated_value
        print(f"- {description}: {calculated_value} points")
        equation_str_parts.append(str(calculated_value))

    # Print the final equation for the top diagnosis
    print(f"\nFinal Equation for the highest scoring diagnosis: {' + '.join(equation_str_parts)} = {final_score}")

diagnose_skin_condition()