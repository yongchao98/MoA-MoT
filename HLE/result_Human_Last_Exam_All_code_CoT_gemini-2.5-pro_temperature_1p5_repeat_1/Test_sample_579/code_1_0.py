def diagnose_skin_condition():
    """
    Analyzes clinical findings to suggest the most likely diagnosis using a scoring system.
    """
    # Define diagnoses and their initial scores
    diagnoses = {
        'A. Malignant Intertrigo': 0,
        'B. Allergic contact dermatitis': 0,
        'C. Hidradenitis Supportiva': 0,
        'D. Atopic dermatitis': 0,
        'E. Psoriasis': 0
    }

    # Clinical findings and risk factors from the case, with associated scores
    # Scores are based on the strength of association with each diagnosis
    evidence = {
        'Multiple Intertriginous Sites': {'C. Hidradenitis Supportiva': 3, 'E. Psoriasis': 2, 'B. Allergic contact dermatitis': 1, 'D. Atopic dermatitis': 1},
        'Purulent Nodules': {'C. Hidradenitis Supportiva': 5},
        'Erythematous Plaques': {'C. Hidradenitis Supportiva': 2, 'E. Psoriasis': 3, 'B. Allergic contact dermatitis': 2, 'D. Atopic dermatitis': 2},
        'Large Bullae': {'C. Hidradenitis Supportiva': 1}, # Can be seen with intense inflammation/abscesses
        'Obesity (BMI > 30)': {'C. Hidradenitis Supportiva': 2, 'E. Psoriasis': 1},
        'Smoking': {'C. Hidradenitis Supportiva': 2, 'E. Psoriasis': 1}
    }

    # Store the equation strings for printing
    equations = {}

    # Calculate scores
    for diagnosis in diagnoses:
        score = 0
        equation_parts = []
        for feature, points_map in evidence.items():
            points = points_map.get(diagnosis, 0)
            score += points
            equation_parts.append(str(points))
        diagnoses[diagnosis] = score
        equations[diagnosis] = " + ".join(equation_parts)

    # Find the best diagnosis
    best_diagnosis = max(diagnoses, key=diagnoses.get)

    # Print the results
    print("Clinical Case Analysis based on a Scoring Model:\n")
    print("Features Scored: [Multiple Sites, Purulent Nodules, Plaques, Bullae, Obesity, Smoking]\n")

    for diagnosis, score in sorted(diagnoses.items(), key=lambda item: item[1], reverse=True):
        print(f"Diagnosis: {diagnosis}")
        print(f"Scoring Equation: {equations[diagnosis]} = {score}")
        print("-" * 30)

    print(f"\nConclusion: The diagnosis with the highest score is '{best_diagnosis}'.")

if __name__ == '__main__':
    diagnose_skin_condition()