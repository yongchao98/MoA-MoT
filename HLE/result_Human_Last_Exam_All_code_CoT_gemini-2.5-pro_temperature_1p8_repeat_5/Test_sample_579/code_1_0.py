import sys
# This is a helper function to avoid breaking the script if run in an environment without a TTY
def tty_print(text):
    if sys.stdout.isatty():
        print(text)

def solve_diagnosis():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    using a programmatic scoring system.
    """
    # Step 1: Define key findings from the patient case.
    patient_findings = {
        "lesions": ["purulent nodules", "large bullae (abscess-like)", "erythematous plaques"],
        "locations": ["axillary folds", "inframammary folds", "inguinal regions"],
        "risk_factors": ["obesity (BMI 39)", "smoking"]
    }

    # Step 2: Define characteristic features for each possible diagnosis.
    diagnostic_criteria = {
        "Malignant Intertrigo": {
            "lesions": ["erythematous plaques", "nodules"],
            "locations": ["axillary folds", "inframammary folds", "inguinal regions"],
            "risk_factors": ["history of carcinoma"], # Patient has ductal carcinoma history.
            "specific_hallmarks": []
        },
        "Allergic contact dermatitis": {
            "lesions": ["erythematous plaques"],
            "locations": ["axillary folds", "inframammary folds", "inguinal regions"],
            "risk_factors": [],
            "specific_hallmarks": []
        },
        "Hidradenitis Suppurativa": {
            "lesions": ["purulent nodules", "erythematous plaques", "large bullae (abscess-like)"],
            "locations": ["axillary folds", "inframammary folds", "inguinal regions"],
            "risk_factors": ["obesity (BMI 39)", "smoking"],
            "specific_hallmarks": ["purulent nodules"]
        },
        "Atopic dermatitis": {
            "lesions": ["erythematous plaques"],
            "locations": ["axillary folds", "inframammary folds", "inguinal regions"], # Flexural involvement is common
            "risk_factors": [],
            "specific_hallmarks": []
        },
        "Psoriasis": {
            "lesions": ["erythematous plaques"], # Characteristic of inverse psoriasis
            "locations": ["axillary folds", "inframammary folds", "inguinal regions"],
            "risk_factors": ["smoking"],
            "specific_hallmarks": []
        }
    }

    scores = {diag: 0 for diag in diagnostic_criteria}
    justification = {diag: [] for diag in diagnostic_criteria}

    tty_print("--- Diagnostic Scoring Analysis ---\n")

    # Step 3: Score each diagnosis based on matching patient findings.
    all_findings = patient_findings["lesions"] + patient_findings["locations"] + patient_findings["risk_factors"]
    
    for finding in all_findings:
        for diagnosis, criteria in diagnostic_criteria.items():
            combined_criteria = criteria["lesions"] + criteria["locations"] + criteria["risk_factors"]
            if finding in combined_criteria:
                scores[diagnosis] += 1
                justification[diagnosis].append(f"+1 for '{finding}'")

    # Step 4: Output the rationale and final equation for each diagnosis.
    for diagnosis, score in scores.items():
        print(f"Diagnosis Evaluation: {diagnosis}")
        equation_parts = justification[diagnosis]
        if not equation_parts:
            print("Score = 0")
        else:
            # We construct a string that looks like an equation: 1 + 1 + ... = score
            equation_str = " + ".join(['1'] * len(equation_parts))
            print(f"Equation: {equation_str} = {score}")
            print("Justification:")
            for reason in equation_parts:
                print(f"  - {reason}")
        print("-" * 20)

    # Step 5: Determine the most likely diagnosis.
    best_diagnosis = max(scores, key=scores.get)
    max_score = scores[best_diagnosis]

    print(f"\n--- Conclusion ---")
    print(f"The diagnosis with the highest score is '{best_diagnosis}' with a score of {max_score}.")
    print("This is because it accounts for the specific lesion type (purulent nodules), the locations (axillary, inframammary, inguinal), and the patient's key risk factors (obesity and smoking).")

solve_diagnosis()
<<<C>>>