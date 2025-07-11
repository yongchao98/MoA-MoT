import sys
import io

# Helper function to print output
def print_analysis(text):
    print(text)

# Step 1: Define the key clinical findings from the case.
patient_profile = {
    "Lesions": {
        "Purulent nodules in inguinal regions": 5,
        "Large bullae in axillary folds": 2,
        "Erythematous plaques in inframammary folds": 1
    },
    "Locations": {
        "Classic intertriginous sites (axilla, inguinal, inframammary)": 3
    },
    "Risk Factors": {
        "Obesity (BMI 39)": 2,
        "Smoking history": 2
    }
}

# Step 2: Define diagnostic criteria for the answer choices.
# We assign weights based on how characteristic a finding is for a disease.
disease_criteria = {
    "A. Malignant Intertrigo": {
        "Erythematous plaques in inframammary folds": 2
    },
    "B. Allergic contact dermatitis": {
        "Erythematous plaques in inframammary folds": 1
    },
    "C. Hidradenitis Suppurativa": {
        "Purulent nodules in inguinal regions": 5, # Hallmark finding
        "Classic intertriginous sites (axilla, inguinal, inframammary)": 3, # Classic distribution
        "Obesity (BMI 39)": 2, # Strong risk factor
        "Smoking history": 2 # Strong risk factor
    },
    "D. Atopic dermatitis": {
        "Erythematous plaques in inframammary folds": 1
    },
    "E. Psoriasis": {
        "Erythematous plaques in inframammary folds": 2
    }
}

print_analysis("### Clinical Case Analysis ###")
print_analysis("\nPatient's Key Findings:")
for category, details in patient_profile.items():
    for finding, _ in details.items():
        print_analysis(f"- {finding}")

print_analysis("\n### Differential Diagnosis Scoring ###")

scores = {}
final_equations = {}

# Step 3: Calculate a score for each potential diagnosis.
for disease, criteria in disease_criteria.items():
    score = 0
    equation_parts = []
    print_analysis(f"\nEvaluating: {disease}")
    
    # Check Lesions
    if "Lesions" in patient_profile:
        for finding, weight in patient_profile["Lesions"].items():
            if finding in criteria:
                score += criteria[finding]
                equation_parts.append(str(criteria[finding]))
                print_analysis(f"  - Finding '{finding}' matches. Adding {criteria[finding]} points.")

    # Check Locations
    if "Locations" in patient_profile:
        for finding, weight in patient_profile["Locations"].items():
            if finding in criteria:
                score += criteria[finding]
                equation_parts.append(str(criteria[finding]))
                print_analysis(f"  - Finding '{finding}' matches. Adding {criteria[finding]} points.")

    # Check Risk Factors
    if "Risk Factors" in patient_profile:
        for finding, weight in patient_profile["Risk Factors"].items():
            if finding in criteria:
                score += criteria[finding]
                equation_parts.append(str(criteria[finding]))
                print_analysis(f"  - Finding '{finding}' matches. Adding {criteria[finding]} points.")

    scores[disease] = score
    final_equations[disease] = " + ".join(equation_parts)
    print_analysis(f"  - Total Score: {score}")

# Step 4: Determine the most likely diagnosis.
most_likely_diagnosis = max(scores, key=scores.get)
winning_score = scores[most_likely_diagnosis]
winning_equation = final_equations[most_likely_diagnosis]

print_analysis("\n### Conclusion ###")
print_analysis(f"The diagnosis with the highest score is '{most_likely_diagnosis}'.")
print_analysis(f"The presence of purulent nodules in characteristic intertriginous sites in a patient with obesity and a smoking history are classic indicators.")
print_analysis("\nFinal Diagnostic Equation:")
# Here we output each number in the final equation as requested.
print(f"{most_likely_diagnosis} Score = {winning_equation} = {winning_score}")
print("<<<C>>>")