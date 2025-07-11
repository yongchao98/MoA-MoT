def diagnose_esophageal_condition():
    """
    Analyzes patient data to determine the most likely esophageal diagnosis
    by assigning and calculating scores for each possibility.
    """
    # Patient findings: 1 for present, 0 for absent, -1 for contradictory
    findings = {
        'Major SCC Risk Factors (Smoking/Alcohol)': 1,
        'Imaging (Wall Thickening/Narrowing)': 1,
        'Normal Endoscopy (No Ulcers/Plaques)': 1,
        'Odynophagia/Chest Pain': 1
    }

    # Scoring matrix based on clinical knowledge
    # [SCC_Risk, Imaging, Normal_Endoscopy, Pain]
    scoring_weights = {
        'A. Streptococcal esophagitis': [0, 0, -2, 1],
        'B. Esophageal adenocarcinoma': [0.5, 1, -2, 1],
        'C. Esophageal squamous cell carcinoma': [3, 2, 1, 1],
        'D. GERD': [0, -2, 0, 1],
        'E. Herpes esophagitis': [0, 0, -2, 1]
    }

    # Calculate final scores
    final_scores = {}
    equations = {}
    for diagnosis, weights in scoring_weights.items():
        score = (weights[0] * findings['Major SCC Risk Factors (Smoking/Alcohol)'] +
                 weights[1] * findings['Imaging (Wall Thickening/Narrowing)'] +
                 weights[2] * findings['Normal Endoscopy (No Ulcers/Plaques)'] +
                 weights[3] * findings['Odynophagia/Chest Pain'])
        final_scores[diagnosis] = score

        # Build the equation string as requested
        # e.g., "Score = (3 * 1) + (2 * 1) + (1 * 1) + (1 * 1) = 7"
        equation = (f"Score = ({weights[0]} * {findings['Major SCC Risk Factors (Smoking/Alcohol)']}) + "
                    f"({weights[1]} * {findings['Imaging (Wall Thickening/Narrowing)']}) + "
                    f"({weights[2]} * {findings['Normal Endoscopy (No Ulcers/Plaques)']}) + "
                    f"({weights[3]} * {findings['Odynophagia/Chest Pain']}) = {score}")
        equations[diagnosis] = equation


    most_likely_diagnosis = max(final_scores, key=final_scores.get)

    print("Diagnostic Score Calculation:")
    print("Based on: [Major SCC Risks, Imaging, Normal Endoscopy, Pain]")
    print("--------------------------------------------------")
    for diagnosis, score in sorted(final_scores.items()):
        print(f"Diagnosis: {diagnosis}")
        print(f"  Calculation: {equations[diagnosis]}")
        print(f"  Final Score: {score}")
        print("")

    print("Conclusion:")
    print(f"The highest score corresponds to the most likely diagnosis.")
    print(f"\nMost Likely Diagnosis: {most_likely_diagnosis}")

diagnose_esophageal_condition()