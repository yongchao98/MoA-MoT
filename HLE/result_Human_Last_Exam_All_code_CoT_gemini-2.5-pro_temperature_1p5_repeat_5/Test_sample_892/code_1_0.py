import sys
import io

# Define a function to encapsulate the logic
def solve_clinical_case():
    """
    This function analyzes a clinical vignette to determine the most likely diagnosis
    by scoring each possibility based on the patient's signs and symptoms.
    """
    # 1. Define patient findings and weights for their significance.
    # The vertebral mass is the most critical finding and receives the highest weight.
    patient_data = {
        "findings": {
            "Chronic Cough": 1,
            "Dyspnea": 1,
            "Acid Reflux": 1,
            "COPD History": 1, # Explains pre-existing conditions but not the new problem
            "Vertebral Mass": 10, # Most critical finding, highly weighted
            "Elevated Creatinine": 2 # Significant finding, could be paraneoplastic
        },
        "present": ["Chronic Cough", "Dyspnea", "Acid Reflux", "COPD History", "Vertebral Mass", "Elevated Creatinine"]
    }

    # 2. Map diagnoses to the findings they typically explain.
    explanation_map = {
        "A. Aspiration pneumonitis": ["Dyspnea", "Acid Reflux"],
        "B. Aspiration pneumonia": ["Chronic Cough", "Dyspnea", "Acid Reflux"],
        "C. Achalasia": ["Chronic Cough", "Acid Reflux"],
        "D. Adenocarcinoma": ["Chronic Cough", "Dyspnea", "Vertebral Mass", "Elevated Creatinine"],
        "E. COPD": ["Chronic Cough", "Dyspnea", "COPD History"]
    }

    # 3. Calculate and print the score for each diagnosis.
    scores = {}
    print("Scoring each diagnosis based on the patient's clinical findings:\n")

    for diagnosis, explained_findings in explanation_map.items():
        score = 0
        equation_parts = []
        for finding in explained_findings:
            if finding in patient_data["present"]:
                weight = patient_data["findings"][finding]
                score += weight
                # Build the equation string part for each finding
                equation_parts.append(f"{finding} ({weight})")
        
        scores[diagnosis] = score
        # Print the final "equation" for the current diagnosis
        final_equation = " + ".join(equation_parts)
        print(f"Diagnosis: {diagnosis}")
        print(f"Score Equation: {final_equation}")
        print(f"Total Score: {score}\n")

    # 4. Determine and print the most likely diagnosis.
    most_likely_diagnosis = max(scores, key=scores.get)

    print("-" * 30)
    print(f"The diagnosis with the highest score is: {most_likely_diagnosis}")
    print("This is because it is the only diagnosis that explains the critical finding of a vertebral mass, suggesting metastatic cancer.")
    print("-" * 30)

# Execute the function
solve_clinical_case()