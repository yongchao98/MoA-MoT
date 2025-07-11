def diagnose_hypoxemia_cause(patient_data, diagnoses):
    """
    Analyzes patient data to determine the most likely cause of hypoxemia
    by scoring potential diagnoses based on clinical fit.
    """
    scores = {}
    print("Analyzing Patient Data:")
    print(f"- Days Post-Op: {patient_data['days_post_op']}")
    print(f"- O2 Saturation: {patient_data['o2_sat']}% on {patient_data['o2_support']}L")
    print(f"- Key Findings: {', '.join(patient_data['findings'])}\n")

    print("Scoring Potential Diagnoses...")
    for key, factors in diagnoses.items():
        # Calculate a simple fit score
        score = factors['timeline_fit'] + factors['symptom_fit'] + factors['context_fit']
        scores[key] = score
        print(f"- Scoring '{key}': Timeline Fit ({factors['timeline_fit']}) + Symptom Fit ({factors['symptom_fit']}) + Context Fit ({factors['context_fit']}) = {score}")

    # Find the diagnosis with the highest score
    most_likely_diagnosis_key = max(scores, key=scores.get)
    most_likely_diagnosis = diagnoses[most_likely_diagnosis_key]
    max_score = scores[most_likely_diagnosis_key]

    print("\n--- Conclusion ---")
    print(f"The most likely diagnosis is: {most_likely_diagnosis_key}")
    print(f"Reasoning: This diagnosis has the highest compatibility score ({max_score}).")
    print(f"It strongly aligns with the timeline of {patient_data['days_post_op']} days post-surgery and explains the severe hypoxemia ({patient_data['o2_sat']}% O2) and bilateral crackles, which are classic signs of sepsis-induced ARDS after major abdominal surgery.")
    
    # Fulfilling the request to show the numbers in a final equation for the top diagnosis
    print("\nFinal Likelihood Score Equation:")
    print(f"Score for '{most_likely_diagnosis_key}' = Timeline Fit ({most_likely_diagnosis['timeline_fit']}) + Symptom Fit ({most_likely_diagnosis['symptom_fit']}) + Context Fit ({most_likely_diagnosis['context_fit']})")


# Patient information from the case
patient_case = {
    "age": 59,
    "days_post_op": 29,
    "o2_sat": 82,
    "o2_support": 3,
    "findings": ["bilateral crackles", "gasping for air (respiratory distress)"]
}

# Scoring matrix for diagnoses (scale 1-5 for fit)
# 1=Poor Fit, 5=Excellent Fit
diagnostic_options = {
    "A. Acute blood transfusion reaction": {"timeline_fit": 1, "symptom_fit": 5, "context_fit": 2},
    "B. Iodine-related reaction":         {"timeline_fit": 1, "symptom_fit": 2, "context_fit": 1},
    "C. Sensitivity reaction":            {"timeline_fit": 2, "symptom_fit": 2, "context_fit": 2},
    "D. Sepsis":                          {"timeline_fit": 5, "symptom_fit": 5, "context_fit": 5},
    "E. Myocyte necrosis":                {"timeline_fit": 2, "symptom_fit": 3, "context_fit": 2},
    "F. Respiratory deconditioning":      {"timeline_fit": 3, "symptom_fit": 1, "context_fit": 2},
    "G. Lung exhaustion":                 {"timeline_fit": 1, "symptom_fit": 1, "context_fit": 1},
    "H. Air pollution sensitivity":       {"timeline_fit": 1, "symptom_fit": 1, "context_fit": 1}
}

# Run the analysis
diagnose_hypoxemia_cause(patient_case, diagnostic_options)