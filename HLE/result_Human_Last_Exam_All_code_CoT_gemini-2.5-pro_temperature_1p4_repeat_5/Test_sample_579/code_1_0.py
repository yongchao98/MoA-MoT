import sys
import io

# Redirect print to a string buffer to control output format
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_diagnosis():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring potential diagnoses against patient findings.
    """
    # Step 1: Define patient's key clinical findings from the vignette
    patient_findings = {
        "Obesity (BMI 39)",
        "Smoker",
        "Type 2 Diabetes Mellitus",
        "Lesions in intertriginous areas (axillary, inframammary, inguinal)",
        "Purulent nodules",
        "Erythematous plaques"
    }

    # Step 2: Define characteristic features for each potential diagnosis
    diagnosis_features = {
        "A. Malignant Intertrigo": {
            "Lesions in intertriginous areas (axillary, inframammary, inguinal)"
            # Note: A vague, non-standard diagnosis, expected to score low.
        },
        "B. Allergic contact dermatitis": {
            "Erythematous plaques",
            "Vesicles",
            "Intense pruritus (itching)"
        },
        "C. Hidradenitis Suppurativa": {
            "Lesions in intertriginous areas (axillary, inframammary, inguinal)",
            "Purulent nodules",
            "Obesity (BMI 39)",  # Strong risk factor
            "Smoker",  # Strong risk factor
            "Type 2 Diabetes Mellitus" # Common comorbidity
        },
        "D. Atopic dermatitis": {
            "Erythematous plaques",
            "Intense pruritus (itching)",
            "Lesions in flexural areas" # Similar to intertriginous, but nodules are not typical
        },
        "E. Psoriasis": {
            "Lesions in intertriginous areas (axillary, inframammary, inguinal)",
            "Erythematous plaques",
            "Obesity (BMI 39)" # Risk factor for inverse psoriasis
        }
    }

    print("--- Diagnostic Scoring Based on Patient Findings ---\n")
    scores = {}

    # Step 3 & 4: Compare patient findings to diagnostic criteria and calculate a score
    for diagnosis, features in diagnosis_features.items():
        matched_features = patient_findings.intersection(features)
        score = len(matched_features)
        scores[diagnosis] = score

        # Print the "equation" for the score calculation for each diagnosis
        print(f"Evaluating: {diagnosis}")
        print(f"  Matching features found: {score}")
        print(f"  Final Score Equation: {score} = { ' + '.join(['1'] * score) if score > 0 else '0' }")
        print("-" * 30)


    # Step 5: Identify the diagnosis with the highest score
    most_likely_diagnosis = max(scores, key=scores.get)

    print("\n--- Conclusion ---")
    print(f"The highest score is {scores[most_likely_diagnosis]} for '{most_likely_diagnosis}'.")
    print("This is based on the strong overlap between the patient's risk factors (obesity, smoking, diabetes), lesion locations (intertriginous folds), and specific morphology (purulent nodules), which are all hallmarks of Hidradenitis Suppurativa.")

# Run the analysis
solve_diagnosis()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())