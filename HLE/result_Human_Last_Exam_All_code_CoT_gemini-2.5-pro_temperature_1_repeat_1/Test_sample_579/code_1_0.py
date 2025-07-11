import collections

def diagnose_skin_condition():
    """
    Analyzes patient findings to suggest the most likely dermatological diagnosis.
    """
    # Key findings from the clinical vignette
    patient_findings = {
        "intertriginous_locations": True, # Axillary, inframammary, inguinal
        "purulent_nodules": True,
        "erythematous_plaques": True,
        "bullae": True,
        "obesity_risk_factor": True, # BMI 39
        "smoking_risk_factor": True
    }

    # Diagnostic criteria and scoring weight for each feature.
    # Positive points for matching features, negative for contradictory ones.
    diagnostic_criteria = {
        "Malignant Intertrigo": {
            "intertriginous_locations": 1, "purulent_nodules": -2, "erythematous_plaques": 1, "bullae": -1
        },
        "Allergic contact dermatitis": {
            "intertriginous_locations": 1, "purulent_nodules": -3, "erythematous_plaques": 1, "bullae": 1
        },
        "Hidradenitis Suppurativa": {
            "intertriginous_locations": 3, "purulent_nodules": 3, "erythematous_plaques": 1, "bullae": 1, # Bullae can be atypical descriptor for abscess
            "obesity_risk_factor": 2, "smoking_risk_factor": 2
        },
        "Atopic dermatitis": {
            "intertriginous_locations": 1, "purulent_nodules": -3, "erythematous_plaques": 1, "bullae": -1
        },
        "Psoriasis": {
            "intertriginous_locations": 2, "purulent_nodules": -3, "erythematous_plaques": 2, "bullae": -2,
            "obesity_risk_factor": 1
        }
    }

    print("Analyzing patient findings against potential diagnoses...\n")

    scores = collections.defaultdict(int)
    equations = collections.defaultdict(list)

    for diagnosis, criteria in diagnostic_criteria.items():
        for feature, present in patient_findings.items():
            if present and feature in criteria:
                score_value = criteria[feature]
                scores[diagnosis] += score_value
                # Build the equation string for each number
                equations[diagnosis].append(f"{feature}: {score_value}")

    # Print the results in a clear format
    print("--- Diagnostic Score Calculation ---")
    for diagnosis, score in sorted(scores.items(), key=lambda item: item[1], reverse=True):
        equation_str = " + ".join(equations[diagnosis]).replace("+ -", "- ")
        print(f"\nDiagnosis: {diagnosis}")
        print(f"Final Score: {score}")
        print(f"Calculation: {equation_str}")


    # Determine the most likely diagnosis
    most_likely_diagnosis = max(scores, key=scores.get)
    print("\n-------------------------------------------")
    print(f"Conclusion: The clinical presentation, including the involvement of multiple intertriginous sites, the presence of purulent nodules, and risk factors like obesity and smoking, most strongly supports the diagnosis of {most_likely_diagnosis}.")
    print("-------------------------------------------")

# Execute the diagnostic function
diagnose_skin_condition()