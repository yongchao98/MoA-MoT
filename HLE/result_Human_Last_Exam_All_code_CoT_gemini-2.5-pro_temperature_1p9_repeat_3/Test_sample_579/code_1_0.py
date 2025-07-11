def diagnose_skin_condition():
    """
    Analyzes patient features to determine the most likely dermatological diagnosis.
    """
    # Patient features from the clinical vignette
    patient_features = {
        "locations": ["axillary folds", "inframammary folds", "inguinal regions"],
        "morphology": ["large bullae", "erythematous plaques", "purulent nodules"],
        "risk_factors": ["obesity (BMI 39)", "smoking"],
        "history": ["ductal carcinoma"]
    }

    # Initialize a score for each potential diagnosis
    scores = {
        "A. Malignant Intertrigo": 0,
        "B. Allergic contact dermatitis": 0,
        "C. Hidradenitis Suppurativa": 0,
        "D. Atopic dermatitis": 0,
        "E. Psoriasis": 0,
    }

    # --- SCORING LOGIC ---

    # 1. Location Analysis: Involvement of multiple classic intertriginous sites
    if len(patient_features["locations"]) == 3:
        # Involvement of axillary, inframammary, and inguinal regions is a hallmark of HS
        scores["C. Hidradenitis Suppurativa"] += 4

    # 2. Morphology Analysis: The type of lesions observed
    if "purulent nodules" in patient_features["morphology"]:
        # Purulent (pus-filled) nodules are highly characteristic of HS
        scores["C. Hidradenitis Suppurativa"] += 5
    if "erythematous plaques" in patient_features["morphology"] and "inframammary folds" in patient_features["locations"]:
        # Plaques in folds can be seen in inverse psoriasis or malignant intertrigo
        scores["E. Psoriasis"] += 1
        # The history of ductal carcinoma makes malignant intertrigo a consideration for this specific site
        if "ductal carcinoma" in patient_features["history"]:
            scores["A. Malignant Intertrigo"] += 2

    # 3. Risk Factor Analysis
    if "obesity (BMI 39)" in patient_features["risk_factors"]:
        # Obesity is a major risk factor for HS
        scores["C. Hidradenitis Suppurativa"] += 2
    if "smoking" in patient_features["risk_factors"]:
        # Smoking is another major risk factor for HS
        scores["C. Hidradenitis Suppurativa"] += 2

    # Determine the diagnosis with the highest score
    most_likely_diagnosis = max(scores, key=scores.get)

    print("--- Diagnostic Analysis ---")
    print("This script evaluates diagnoses based on a scoring system.")
    print("\nLikelihood scores based on patient features:")
    # Using a loop to print each score, fulfilling the requirement to output each number
    for diagnosis, score in scores.items():
        print(f"Final Score for {diagnosis}: {score}")

    print(f"\nConclusion: The clinical presentation strongly suggests the diagnosis with the highest score.")
    print(f"Most Likely Diagnosis: {most_likely_diagnosis}")


# Execute the diagnostic function
diagnose_skin_condition()