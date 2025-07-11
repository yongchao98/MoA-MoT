def diagnose_skin_condition():
    """
    This function analyzes a clinical case by scoring potential diagnoses
    based on the patient's signs, symptoms, and risk factors.
    """

    # Initialize scores for each potential diagnosis
    scores = {
        "A. Malignant Intertrigo": 0,
        "B. Allergic contact dermatitis": 0,
        "C. Hidradenitis Supportiva": 0,
        "D. Atopic dermatitis": 0,
        "E. Psoriasis": 0,
    }

    print("Step 1: Analyzing the primary lesion type - Purulent Nodules.")
    # Purulent (pus-filled) nodules are a hallmark of Hidradenitis Suppurativa (HS).
    # This is the most specific finding.
    print("Finding: Purulent nodules in the inguinal regions.")
    print("Rationale: This is a classic and defining feature of HS.")
    scores["C. Hidradenitis Supportiva"] += 3
    print("Current Scores Equation:")
    print(f"C. Hidradenitis Supportiva = {scores['C. Hidradenitis Supportiva']}")

    print("\nStep 2: Analyzing the location of the lesions.")
    # The locations are classic for diseases affecting skin folds.
    print("Finding: Lesions in axillary, inframammary, and inguinal folds.")
    print("Rationale: This distribution is highly characteristic of HS. It is also seen in Inverse Psoriasis.")
    scores["C. Hidradenitis Supportiva"] += 2
    scores["E. Psoriasis"] += 1
    print("Current Scores Equation:")
    print(f"C. Hidradenitis Supportiva = {scores['C. Hidradenitis Supportiva'] - 2} + 2 = {scores['C. Hidradenitis Supportiva']}")
    print(f"E. Psoriasis = {scores['E. Psoriasis'] - 1} + 1 = {scores['E. Psoriasis']}")


    print("\nStep 3: Analyzing patient risk factors.")
    # Obesity (BMI 39) and smoking are major risk factors for HS.
    print("Finding: Patient has a BMI of 39 (obesity) and is a smoker.")
    print("Rationale: Obesity and smoking are two of the strongest risk factors for developing HS.")
    scores["C. Hidradenitis Supportiva"] += 2
    print("Current Scores Equation:")
    print(f"C. Hidradenitis Supportiva = {scores['C. Hidradenitis Supportiva'] - 2} + 2 = {scores['C. Hidradenitis Supportiva']}")

    print("\nStep 4: Analyzing other lesion types.")
    # Erythematous plaques and bullae can be seen in several conditions, but fit the inflammatory nature of HS.
    # Large, painful abscesses in HS can be described as bullae.
    print("Finding: Large bullae and erythematous plaques.")
    print("Rationale: These are consistent with the inflammation and abscess formation in HS.")
    scores["C. Hidradenitis Supportiva"] += 1
    scores["E. Psoriasis"] += 1 # Plaques are also a feature of psoriasis
    print("Current Scores Equation:")
    print(f"C. Hidradenitis Supportiva = {scores['C. Hidradenitis Supportiva'] - 1} + 1 = {scores['C. Hidradenitis Supportiva']}")
    print(f"E. Psoriasis = {scores['E. Psoriasis'] - 1} + 1 = {scores['E. Psoriasis']}")

    print("\n--- Final Diagnostic Score Calculation ---")
    for diagnosis, score in scores.items():
        print(f"Total score for {diagnosis}: {score}")

    # Determine the most likely diagnosis
    best_diagnosis = max(scores, key=scores.get)

    print(f"\nConclusion: The diagnosis with the highest score is '{best_diagnosis}'.")

if __name__ == "__main__":
    diagnose_skin_condition()