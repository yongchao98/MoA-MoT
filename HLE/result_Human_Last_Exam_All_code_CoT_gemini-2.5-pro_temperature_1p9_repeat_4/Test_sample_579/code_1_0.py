def solve_medical_case():
    """
    This function analyzes a clinical case by scoring potential diagnoses
    based on the patient's symptoms, history, and risk factors.
    """

    # --- Patient Data Points from the Case ---
    findings = {
        'age': 64,
        'bmi': 39,
        'is_smoker': True,
        'history_of_cancer': True,
        'locations': ['axillary folds', 'inframammary folds', 'inguinal regions'],
        'lesions': ['large bullae', 'erythematous plaques', 'purulent nodules']
    }

    # --- Initialize scores for each diagnosis ---
    scores = {
        'A. Malignant Intertrigo': 0,
        'B. Allergic contact dermatitis': 0,
        'C. Hidradenitis Supportiva': 0,
        'D. Atopic dermatitis': 0,
        'E. Psoriasis': 0
    }

    print("Analyzing clinical findings to determine the most likely diagnosis...")
    print("-" * 30)

    # --- Scoring Logic ---

    # 1. Evaluate hallmark lesion: Purulent nodules
    if 'purulent nodules' in findings['lesions']:
        scores['C. Hidradenitis Supportiva'] += 3
        print("Finding: 'Purulent nodules' strongly suggests Hidradenitis Supportiva. Score +3")

    # 2. Evaluate locations: Classic intertriginous sites
    classic_hs_sites = ['axillary', 'inguinal', 'inframammary']
    affected_hs_sites = sum(1 for site in findings['locations'] if any(hs_site in site for hs_site in classic_hs_sites))
    if affected_hs_sites >= 2:
        scores['C. Hidradenitis Supportiva'] += 3
        print(f"Finding: Lesions in {affected_hs_sites} classic sites (axilla, inguinal, inframammary) strongly suggests Hidradenitis Supportiva. Score +3")

    # 3. Evaluate risk factors: Obesity (BMI > 30)
    if findings['bmi'] > 30:
        scores['C. Hidradenitis Supportiva'] += 2
        scores['E. Psoriasis'] += 1
        print(f"Finding: High BMI ({findings['bmi']}) is a major risk factor for HS (Score +2) and a minor one for Psoriasis (Score +1).")

    # 4. Evaluate risk factors: Smoking
    if findings['is_smoker']:
        scores['C. Hidradenitis Supportiva'] += 2
        print("Finding: Smoking is a major risk factor for Hidradenitis Supportiva. Score +2")

    # 5. Evaluate history: Cancer
    if findings['history_of_cancer']:
        scores['A. Malignant Intertrigo'] += 2
        print("Finding: History of cancer is a risk factor for Malignant Intertrigo. Score +2")
    
    # 6. Evaluate non-specific lesions: Plaques
    if 'erythematous plaques' in findings['lesions']:
        scores['A. Malignant Intertrigo'] += 1
        scores['E. Psoriasis'] += 1
        scores['D. Atopic dermatitis'] += 1
        scores['B. Allergic contact dermatitis'] += 1
        print("Finding: 'Erythematous plaques' is a non-specific finding, adding minor scores to several possibilities. Score +1 each")


    print("-" * 30)
    print("Final Scores:")
    for diagnosis, score in scores.items():
        print(f"{diagnosis}: {score}")
    print("-" * 30)

    # Determine the most likely diagnosis
    best_diagnosis = max(scores, key=scores.get)
    print(f"The highest scoring diagnosis is: {best_diagnosis}")

if __name__ == "__main__":
    solve_medical_case()
    print("\nThe correct answer choice is C, which corresponds to Hidradenitis Supportiva.")
    print('<<<C>>>')