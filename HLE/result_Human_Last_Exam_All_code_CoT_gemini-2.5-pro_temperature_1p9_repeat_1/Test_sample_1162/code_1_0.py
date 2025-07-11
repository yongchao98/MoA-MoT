import operator

def solve_diagnosis():
    """
    This script analyzes patient symptoms to determine the most likely diagnosis
    from a list of options by using a scoring system.
    """
    
    # Define the key symptoms from the clinical vignette
    patient_symptoms = {
        "age_2_years": True,
        "hypertension": True,
        "growth_delay": True,
        "aniridia": True,       # A very specific and high-yield finding
        "developmental_delay": True,
        "pelvic_mass": True
    }
    
    # Scoring matrix: points assigned for each symptom corresponding to a diagnosis.
    # The aniridia-nephroblastoma link is a classic association (WAGR syndrome)
    # and is given a very high score to reflect its diagnostic importance.
    scoring_matrix = {
        "A. Germ cell tumor":  {"pelvic_mass": 2, "age_2_years": 1},
        "B. Astrocytoma":      {},
        "C. Neuroblastoma":    {"age_2_years": 3, "hypertension": 3, "pelvic_mass": 3, "growth_delay": 1},
        "D. Nephroblastoma":   {"age_2_years": 3, "hypertension": 3, "pelvic_mass": 3, "aniridia": 10, "developmental_delay": 2, "growth_delay": 1},
        "E. Ewing sarcoma":    {"pelvic_mass": 2}
    }
    
    print("Calculating diagnostic scores based on patient's symptoms:")
    print("-" * 60)
    
    final_scores = {}
    
    # Calculate the total score for each diagnosis
    for diagnosis, points_map in scoring_matrix.items():
        total_score = 0
        equation_parts = []
        for symptom, present in patient_symptoms.items():
            if present and symptom in points_map:
                score = points_map[symptom]
                total_score += score
                equation_parts.append(f"{score} (from {symptom})")
        
        final_scores[diagnosis] = total_score
        
        # Print the breakdown of the calculation for each diagnosis
        if not equation_parts:
            equation_str = "0"
        else:
            equation_str = " + ".join(equation_parts)
            
        print(f"Diagnosis: {diagnosis}")
        print(f"Score calculation: {equation_str} = {total_score}\n")

    # Determine the most likely diagnosis
    if not final_scores:
        most_likely_diagnosis = "Insufficient data"
    else:
        most_likely_diagnosis = max(final_scores.items(), key=operator.itemgetter(1))[0]
    
    print("-" * 60)
    print(f"Conclusion: The diagnosis with the highest score is '{most_likely_diagnosis}'.")
    print("The combination of aniridia, a pelvic (renal) mass, and developmental delay in a young child is classic for WAGR syndrome, which is strongly associated with Nephroblastoma (Wilms tumor).")

solve_diagnosis()