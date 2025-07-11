def diagnose_pediatric_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring each option against the patient's key findings.
    """
    # Step 1: Define the key findings from the clinical vignette.
    patient_findings = {
        "Age 2 years",
        "Hypertension",
        "Aniridia",
        "Developmental Delay",
        "Pelvic/Abdominal Mass"
    }

    # Step 2: Define the characteristic features for each possible diagnosis.
    # This represents a simplified medical knowledge base.
    diagnostic_features = {
        "A. Germ cell tumor": {"Pelvic/Abdominal Mass", "Age 2 years"},
        "B. Astrocytoma": {"Age 2 years"}, # Brain tumor, doesn't match mass location
        "C. Neuroblastoma": {"Age 2 years", "Hypertension", "Pelvic/Abdominal Mass"},
        "D. Nephroblastoma": {"Age 2 years", "Hypertension", "Aniridia", "Developmental Delay", "Pelvic/Abdominal Mass"},
        "E. Ewing sarcoma": {"Pelvic/Abdominal Mass"}
    }

    print("Analyzing patient findings against potential diagnoses...\n")

    best_match = None
    max_score = -1
    score_details = {}

    # Step 3: Score each diagnosis based on matching features.
    for diagnosis, features in diagnostic_features.items():
        matched_items = patient_findings.intersection(features)
        score = len(matched_items)
        score_details[diagnosis] = {"score": score, "matched_items": matched_items}

        if score > max_score:
            max_score = score
            best_match = diagnosis

    # Step 4: Print the results, including the reasoning for the best match.
    print("--- Diagnostic Scoring ---")
    for diagnosis, details in score_details.items():
        print(f"Diagnosis: {diagnosis:<25} | Score: {details['score']}")
    print("--------------------------\n")
    
    print(f"The most likely diagnosis is: {best_match}")
    print("This is because it has the highest number of matching features with the patient's presentation.")
    
    # Step 5: Fulfill the request to show the final equation with numbers.
    # Each matched feature contributes 1 point to the score.
    best_match_info = score_details[best_match]
    matched_list = sorted(list(best_match_info["matched_items"]))
    
    # Create the equation string: 1 + 1 + 1... for each matched feature.
    equation_numbers = ["1"] * len(matched_list)
    equation_str = " + ".join(equation_numbers)
    
    print("\n--- Final Score Calculation for Best Match ---")
    print(f"Matching Features: {', '.join(matched_list)}")
    print(f"Equation: {best_match.split('.')[1].strip()} Score = {equation_str}")
    print(f"Final Score = {best_match_info['score']}")


if __name__ == "__main__":
    diagnose_pediatric_case()