def solve_diagnosis():
    """
    This script analyzes clinical findings to determine the most likely diagnosis
    by using a weighted scoring system.
    """

    # --- Patient's Clinical Findings ---
    # We assign a point value to each finding based on its diagnostic significance.
    # Aniridia is given a high score because of its strong association with a specific diagnosis.
    findings = {
        'Aniridia': 5,
        'Pelvic Mass': 3,
        'Hypertension': 2,
        'Developmental Delay': 2,
        'Age 2yo': 1
    }
    
    # --- Potential Diagnoses and their associated symptoms ---
    # We define which findings are typically associated with each disease.
    associations = {
        'A. Germ cell tumor': ['Pelvic Mass', 'Age 2yo'],
        'B. Astrocytoma': [], # Brain tumor, does not fit the pelvic mass
        'C. Neuroblastoma': ['Pelvic Mass', 'Hypertension', 'Age 2yo'],
        'D. Nephroblastoma': ['Aniridia', 'Pelvic Mass', 'Hypertension', 'Developmental Delay', 'Age 2yo'],
        'E. Ewing sarcoma': ['Pelvic Mass', 'Age 2yo']
    }

    print("Evaluating diagnostic possibilities based on a scoring system...\n")
    
    final_scores = {}
    best_diagnosis = ''
    max_score = -1

    # --- Calculate and Print Scores ---
    # The script iterates through each diagnosis, calculates a score, and prints the equation.
    for diagnosis, associated_findings in associations.items():
        score = 0
        equation_parts = []
        for finding in associated_findings:
            point_value = findings.get(finding, 0)
            score += point_value
            equation_parts.append(str(point_value))
        
        final_scores[diagnosis] = score
        
        # We format and print the equation for each diagnosis as requested.
        equation = " + ".join(equation_parts) if equation_parts else "0"
        print(f"Diagnosis: {diagnosis}")
        print(f"Equation based on associated findings: {equation} = {score}")
        print("-" * 20)

        if score > max_score:
            max_score = score
            best_diagnosis = diagnosis

    print(f"\nConclusion: The highest scoring diagnosis is {best_diagnosis} with {max_score} points.")
    print("The combination of aniridia, a pelvic mass (Wilms tumor), developmental delay, and hypertension strongly points to Nephroblastoma as part of WAGR syndrome.")
    
    # Extract the letter for the final answer format
    final_answer_letter = best_diagnosis.split('.')[0]
    print(f"<<<{final_answer_letter}>>>")

solve_diagnosis()