def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring each option based on key patient findings.
    """
    # Key findings from the patient's case
    patient_findings = {
        "pelvic_mass",
        "hypertension",
        "aniridia",
        "developmental_delay"  # Represented by delayed speech
    }

    # Associated findings for each diagnosis choice
    disease_profiles = {
        "A. Germ cell tumor": {"pelvic_mass"},
        "B. Astrocytoma": {"developmental_delay"},
        "C. Neuroblastoma": {"pelvic_mass", "hypertension"},
        "D. Nephroblastoma": {"pelvic_mass", "hypertension", "aniridia", "developmental_delay"},
        "E. Ewing sarcoma": {"mass"} # Generic mass, less specific
    }

    print("Analyzing patient findings against possible diagnoses:\n")

    diagnosis_scores = {}
    best_diagnosis = ""
    max_score = -1

    # Calculate a score for each diagnosis based on matching findings
    for diagnosis, profile in disease_profiles.items():
        # Find the intersection of patient findings and disease profile
        matching_findings = patient_findings.intersection(profile)
        score = len(matching_findings)
        diagnosis_scores[diagnosis] = score

        # Build the equation string for printing
        # e.g., "Score = 1 (pelvic_mass) + 1 (hypertension) + 1 (aniridia) + 1 (developmental_delay)"
        equation_parts = [f"1 ({finding})" for finding in sorted(list(matching_findings))]
        if not equation_parts:
            equation_string = "Score = 0"
        else:
            equation_string = "Score = " + " + ".join(equation_parts)
        
        print(f"Diagnosis: {diagnosis}")
        print(f"Matching Findings: {sorted(list(matching_findings)) if matching_findings else 'None'}")
        print(f"{equation_string} = {score}\n")

        if score > max_score:
            max_score = score
            best_diagnosis = diagnosis

    print("-" * 30)
    print(f"Conclusion: The diagnosis with the highest score is '{best_diagnosis}' with a score of {max_score}.")
    print("The patient's presentation with a pelvic mass, aniridia, developmental delay (speech), and hypertension is classic for Nephroblastoma (Wilms tumor) in the context of WAGR syndrome.")

# Run the analysis
solve_clinical_case()