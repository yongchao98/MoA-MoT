import sys
import io

# Redirect stdout to capture the print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def run_diagnostic_analysis():
    """
    Analyzes patient findings to suggest a diagnosis using a scoring model.
    """
    # Key findings from the clinical case
    patient_findings = {
        "locations": ["axillary", "inframammary", "inguinal"],
        "lesions": ["bullae", "plaques", "purulent nodules"],
        "risk_factors": ["obesity", "smoking", "diabetes"]
    }

    # Scoring criteria for each diagnosis based on clinical knowledge
    # Positive points for matching features, negative for contradictory ones
    diagnosis_criteria = {
        "A. Malignant Intertrigo": {
            "locations": {"inframammary": 1},
            "lesions": {"plaques": 1, "purulent nodules": -2},
            "risk_factors": {}
        },
        "B. Allergic contact dermatitis": {
            "locations": {"axillary": 1, "inframammary": 1, "inguinal": 1},
            "lesions": {"plaques": 1, "purulent nodules": -3},
            "risk_factors": {}
        },
        "C. Hidradenitis Supportiva": {
            "locations": {"axillary": 3, "inframammary": 3, "inguinal": 3},
            "lesions": {"purulent nodules": 5, "plaques": 2, "bullae": 1},
            "risk_factors": {"obesity": 2, "smoking": 2, "diabetes": 1}
        },
        "D. Atopic dermatitis": {
            "locations": {"axillary": 1, "inframammary": 1},
            "lesions": {"plaques": 1, "purulent nodules": -3},
            "risk_factors": {}
        },
        "E. Psoriasis": {
            "locations": {"axillary": 2, "inframammary": 2, "inguinal": 2},
            "lesions": {"plaques": 3, "purulent nodules": -4},
            "risk_factors": {"obesity": 1, "smoking": 1}
        }
    }

    # Initialize scores and calculation breakdown
    scores = {diag: 0 for diag in diagnosis_criteria}
    score_calculation_parts = {diag: [] for diag in diagnosis_criteria}

    # --- Scoring Logic ---
    # Score based on locations
    for finding in patient_findings["locations"]:
        for diag, criteria in diagnosis_criteria.items():
            if finding in criteria.get("locations", {}):
                score_value = criteria["locations"][finding]
                scores[diag] += score_value
                score_calculation_parts[diag].append(str(score_value))

    # Score based on lesions
    for finding in patient_findings["lesions"]:
        for diag, criteria in diagnosis_criteria.items():
            if finding in criteria.get("lesions", {}):
                score_value = criteria["lesions"][finding]
                scores[diag] += score_value
                score_calculation_parts[diag].append(str(score_value))

    # Score based on risk factors
    for finding in patient_findings["risk_factors"]:
        for diag, criteria in diagnosis_criteria.items():
            if finding in criteria.get("risk_factors", {}):
                score_value = criteria["risk_factors"][finding]
                scores[diag] += score_value
                score_calculation_parts[diag].append(str(score_value))

    # --- Output Results ---
    print("Clinical Feature Scoring for Differential Diagnosis:\n")
    for diag, score in sorted(scores.items(), key=lambda item: item[1], reverse=True):
        # Format the equation string, handling negative numbers
        equation_str = " + ".join(score_calculation_parts[diag]).replace("+ -", "- ")
        print(f"{diag}")
        print(f"Score Equation: {equation_str} = {score}")
        print("-" * 30)

    # Determine the most likely diagnosis
    best_diagnosis = max(scores, key=scores.get)
    print(f"\nConclusion: The highest score indicates the most likely diagnosis is {best_diagnosis}.")

# Run the analysis
run_diagnostic_analysis()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
