def diagnose_patient():
    """
    Analyzes patient symptoms to determine the most likely diagnosis
    using a weighted scoring system.
    """
    # Patient's key clinical findings
    symptoms = {
        "Pelvic Mass": 1,
        "Hypertension": 1,
        "Aniridia": 3,  # This is a highly specific finding, so it gets a higher weight.
        "Developmental Delay": 1,
        "Age_2_years": 1
    }

    # Potential diagnoses and their association with the findings
    # A value of 1 means associated, 0 means not typically associated.
    associations = {
        "Germ cell tumor": {"Pelvic Mass": 1, "Hypertension": 0, "Aniridia": 0, "Developmental Delay": 0, "Age_2_years": 1},
        "Astrocytoma": {"Pelvic Mass": 0, "Hypertension": 0, "Aniridia": 0, "Developmental Delay": 0, "Age_2_years": 0},
        "Neuroblastoma": {"Pelvic Mass": 1, "Hypertension": 1, "Aniridia": 0, "Developmental Delay": 0, "Age_2_years": 1},
        "Nephroblastoma": {"Pelvic Mass": 1, "Hypertension": 1, "Aniridia": 1, "Developmental Delay": 1, "Age_2_years": 1},
        "Ewing sarcoma": {"Pelvic Mass": 1, "Hypertension": 0, "Aniridia": 0, "Developmental Delay": 0, "Age_2_years": 0}
    }

    # Calculate scores for each diagnosis
    scores = {}
    for diagnosis, assoc_symptoms in associations.items():
        total_score = 0
        for symptom, weight in symptoms.items():
            if assoc_symptoms.get(symptom, 0) == 1:
                total_score += weight
        scores[diagnosis] = total_score

    # Find the diagnosis with the highest score
    most_likely_diagnosis = max(scores, key=scores.get)
    
    print("Evaluating diagnoses based on patient's symptoms...\n")
    
    # Print the "equation" for the most likely diagnosis
    print(f"Scoring for {most_likely_diagnosis}:")
    equation_parts = []
    final_score = 0
    for symptom, weight in symptoms.items():
        if associations[most_likely_diagnosis].get(symptom, 0) == 1:
            points = weight
            print(f"- Symptom '{symptom}' adds {points} point(s).")
            equation_parts.append(str(points))
            final_score += points
    
    equation_str = " + ".join(equation_parts)
    print(f"\nFinal Equation: {equation_str} = {final_score}")

    print("\n---")
    print(f"Conclusion: The highest score is for {most_likely_diagnosis}.")
    print("This is because the combination of an abdominal/pelvic mass, hypertension, aniridia, and developmental delay is characteristic of WAGR syndrome, which is strongly associated with Nephroblastoma (Wilms tumor).")

diagnose_patient()