import sys
import io

# Redirect stdout to capture the print output for the final answer
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def evaluate_hypoxemia_case():
    """
    Analyzes a clinical case of post-operative hypoxemia to determine the most likely cause.
    This script scores potential diagnoses based on key features of the case.
    """

    # --- Case Data ---
    # Timeline: 29 days post-op is late for acute reactions, prime time for infections.
    # Clinical Picture: Severe hypoxemia + bilateral crackles = classic ARDS.
    # Risk Factors: Whipple procedure has a high risk of intra-abdominal infection.
    case_findings = {
        "timeline_fit": {
            "Sepsis": 3, "Acute blood transfusion reaction": 0, "Iodine-related reaction": 0,
            "Sensitivity reaction": 1, "Myocyte necrosis": 1, "Respiratory deconditioning": 2,
            "Lung exhaustion": 1, "Air pollution sensitivity": 0
        },
        "symptoms_fit": {
            "Sepsis": 3, "Acute blood transfusion reaction": 1, "Iodine-related reaction": 1,
            "Sensitivity reaction": 0, "Myocyte necrosis": 2, "Respiratory deconditioning": 0,
            "Lung exhaustion": 0, "Air pollution sensitivity": 0
        },
        "risk_factor_fit": {
            "Sepsis": 3, "Acute blood transfusion reaction": 1, "Iodine-related reaction": 0,
            "Sensitivity reaction": 1, "Myocyte necrosis": 1, "Respiratory deconditioning": 2,
            "Lung exhaustion": 1, "Air pollution sensitivity": 0
        }
    }

    diagnoses = {
        "A": "Acute blood transfusion reaction",
        "B": "Iodine-related reaction",
        "C": "Sensitivity reaction",
        "D": "Sepsis",
        "E": "Myocyte necrosis",
        "F": "Respiratory deconditioning",
        "G": "Lung exhaustion",
        "H": "Air pollution sensitivity"
    }

    print("Analyzing potential causes based on a scoring system (0=Poor Fit, 3=Excellent Fit)...\n")

    results = {}
    for key, diagnosis_name in diagnoses.items():
        # Retrieve scores for the current diagnosis
        score_timeline = case_findings["timeline_fit"].get(diagnosis_name, 0)
        score_symptoms = case_findings["symptoms_fit"].get(diagnosis_name, 0)
        score_risk = case_findings["risk_factor_fit"].get(diagnosis_name, 0)

        # Calculate the total score
        total_score = score_timeline + score_symptoms + score_risk

        # Store the result
        results[key] = total_score

        # Print the scoring equation as requested
        print(f"Diagnosis ({key}) {diagnosis_name}:")
        print(f"Score Equation: {score_timeline} (Timeline) + {score_symptoms} (Symptoms) + {score_risk} (Risk Factors) = {total_score}\n")

    # Find the diagnosis with the highest score
    most_likely_key = max(results, key=results.get)
    most_likely_diagnosis = diagnoses[most_likely_key]
    highest_score = results[most_likely_key]

    print("-" * 40)
    print(f"Conclusion: The highest scoring diagnosis is '{most_likely_diagnosis}' with a score of {highest_score}.")
    print("This is because the patient's presentation of ARDS at 29 days post-Whipple procedure is classic for Sepsis originating from a post-operative complication like an abscess.")
    print("-" * 40)

    # This will be captured and printed at the end
    print(f"<<<{most_likely_key}>>>")

# Execute the function
evaluate_hypoxemia_case()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())