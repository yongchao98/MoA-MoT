def solve_diagnosis():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    using a scoring system.
    """
    # Key findings from the clinical vignette
    patient_findings = {
        "Age (2 years)": 1,
        "Hypertension": 1,
        "Aniridia": 1,
        "Abdominal/Pelvic Mass": 1,
        "Developmental Delay": 1,
    }

    # How each diagnosis typically aligns with the key findings (1 for match, 0 for no match)
    # This data is based on established medical knowledge.
    diagnoses_info = {
        "A. Germ cell tumor":       {"Age (2 years)": 1, "Hypertension": 0, "Aniridia": 0, "Abdominal/Pelvic Mass": 1, "Developmental Delay": 0},
        "B. Astrocytoma":           {"Age (2 years)": 0, "Hypertension": 0, "Aniridia": 0, "Abdominal/Pelvic Mass": 0, "Developmental Delay": 0},
        "C. Neuroblastoma":         {"Age (2 years)": 1, "Hypertension": 1, "Aniridia": 0, "Abdominal/Pelvic Mass": 1, "Developmental Delay": 0},
        "D. Nephroblastoma":        {"Age (2 years)": 1, "Hypertension": 1, "Aniridia": 1, "Abdominal/Pelvic Mass": 1, "Developmental Delay": 1},
        "E. Ewing sarcoma":         {"Age (2 years)": 0, "Hypertension": 0, "Aniridia": 0, "Abdominal/Pelvic Mass": 1, "Developmental Delay": 0},
    }

    print("Evaluating potential diagnoses based on a scoring system.")
    print("Each diagnosis is scored based on how many key patient findings it explains.\n")
    print(f"Patient's Key Findings: {', '.join(patient_findings.keys())}\n")

    best_diagnosis = ""
    max_score = -1

    # This loop generates the "equation" of scores for each diagnosis
    for diagnosis, features in diagnoses_info.items():
        score = 0
        equation_parts = []
        for finding_key in patient_findings.keys():
            # Get the score for this feature (0 if not present in the diagnosis info)
            match_score = features.get(finding_key, 0)
            score += match_score
            equation_parts.append(str(match_score))

        # This print statement fulfills the requirement to "output each number in the final equation"
        equation_str = " + ".join(equation_parts)
        print(f"Analysis for {diagnosis}:")
        print(f"  Score Calculation (for each key finding): {equation_str} = {score}")

        if score > max_score:
            max_score = score
            best_diagnosis = diagnosis
        print("-" * 30)

    print("\n--- Conclusion ---")
    print(f"The diagnosis with the highest score ({max_score}) is: {best_diagnosis}")
    print("\nReasoning: The patient's constellation of symptoms, especially the presence of aniridia with an abdominal mass, is classic for WAGR syndrome, which involves Wilms tumor (Nephroblastoma).")

solve_diagnosis()