def diagnose_patient():
    """
    Analyzes clinical findings to determine the most likely diagnosis
    by scoring each option based on symptom correlation.
    """
    # Patient's clinical findings from the case description
    patient_symptoms = {
        "pelvic_mass": 1,
        "hypertension": 1,
        "aniridia": 5,  # Given a high weight due to its high specificity for WAGR syndrome
        "developmental_delay": 2, # Also a key feature of WAGR
        "anemia": 1,
        "failure_to_thrive": 1,
        "typical_age_group": 1 # 2 years old is a typical age for several pediatric tumors
    }

    # Symptom profiles for each potential diagnosis
    diagnosis_profiles = {
        "A. Germ cell tumor": {"pelvic_mass": 1, "typical_age_group": 1},
        "B. Astrocytoma": {"developmental_delay": 1}, # Brain tumor, does not cause pelvic mass
        "C. Neuroblastoma": {"pelvic_mass": 1, "hypertension": 1, "failure_to_thrive": 1, "typical_age_group": 1},
        "D. Nephroblastoma": {"pelvic_mass": 1, "hypertension": 1, "aniridia": 1, "developmental_delay": 1, "anemia": 1, "failure_to_thrive": 1, "typical_age_group": 1},
        "E. Ewing sarcoma": {"pelvic_mass": 1}
    }

    # Calculate and store the score for each diagnosis
    scores = {}
    print("--- Diagnostic Analysis ---")
    for diagnosis, profile in diagnosis_profiles.items():
        score = 0
        # The 'equation' is a sum of weights for matching symptoms
        equation_str_parts = []
        for symptom, weight in patient_symptoms.items():
            if symptom in profile:
                score += weight
                equation_str_parts.append(str(weight))
        
        scores[diagnosis] = score
        # The prompt asks to output each number in the final equation
        equation_str = " + ".join(equation_str_parts) if equation_str_parts else "0"
        print(f"Score for {diagnosis}: {score} (Calculation: {equation_str})")

    # Determine the best diagnosis based on the highest score
    best_diagnosis = max(scores, key=scores.get)
    print("\n--- Conclusion ---")
    print(f"The highest score belongs to: {best_diagnosis}")
    print("This is because Nephroblastoma (Wilms' tumor) is a key component of WAGR syndrome, which perfectly explains the unique combination of a pelvic mass, hypertension, aniridia, and developmental delay.")

if __name__ == "__main__":
    diagnose_patient()