def diagnose_esophageal_condition():
    """
    This script analyzes a clinical vignette to determine the most likely diagnosis.
    It uses a scoring system based on the patient's risk factors, symptoms, and test results.
    """

    # Patient findings from the case
    findings = {
        "Heavy Smoking (>20 pack-years)": {"SCC": 3, "Adeno": 1, "GERD": 0, "Strep": 0, "Herpes": 0},
        "Alcohol Use Disorder": {"SCC": 3, "Adeno": 0, "GERD": 0, "Strep": 0, "Herpes": 0},
        "Severe Pain/Odynophagia": {"SCC": 2, "Adeno": 2, "GERD": 1, "Strep": 2, "Herpes": 2},
        "Imaging (Wall Thickening/Narrowing)": {"SCC": 2, "Adeno": 2, "GERD": 0, "Strep": 1, "Herpes": 1},
        "Labs (Inflammation)": {"SCC": 1, "Adeno": 1, "GERD": 0, "Strep": 2, "Herpes": 1},
        "Endoscopy (No Ulcers/Plaques/Erythema)": {"SCC": 1, "Adeno": -2, "GERD": -3, "Strep": -3, "Herpes": -3}
    }

    # Rationale for the most tricky finding
    print("Clinical Analysis:\n")
    print("The key to this diagnosis is interpreting the 'normal' endoscopy in the context of the other findings.")
    print("While most cancers are visible, SCC can grow underneath the surface (infiltratively),")
    print("causing wall thickening without a visible ulcer or mass initially.\n")
    
    # Initialize scores
    scores = {
        "A. Streptococcal esophagitis": 0,
        "B. Esophageal adenocarcinoma": 0,
        "C. Esophageal squamous cell carcinoma": 0,
        "D. GERD": 0,
        "E. Herpes esophagitis": 0
    }
    
    # Abbreviation mapping
    abbr_map = {
        "C. Esophageal squamous cell carcinoma": "SCC",
        "B. Esophageal adenocarcinoma": "Adeno",
        "D. GERD": "GERD",
        "A. Streptococcal esophagitis": "Strep",
        "E. Herpes esophagitis": "Herpes"
    }
    
    # Calculate and print scores for each diagnosis
    print("Scoring each diagnosis based on the patient's presentation:\n")
    for diagnosis, score in scores.items():
        print(f"--- {diagnosis} ---")
        total_score = 0
        equation_parts = []
        for finding_name, points_map in findings.items():
            abbr = abbr_map[diagnosis]
            point = points_map[abbr]
            total_score += point
            equation_parts.append(str(point))
        
        scores[diagnosis] = total_score
        equation_str = " + ".join(equation_parts).replace("+ -", "- ")
        print(f"Calculation: {equation_str} = {total_score}")
        print(f"Final Score: {total_score}\n")

    # Determine the most likely diagnosis
    most_likely_diagnosis = max(scores, key=scores.get)
    
    print("--- Conclusion ---")
    print(f"The diagnosis with the highest score is: {most_likely_diagnosis}")
    print("The patient's significant risk factors (smoking, alcohol) and clinical findings (infiltrative growth on imaging with normal-appearing mucosa) strongly point to this diagnosis.")

diagnose_esophageal_condition()
<<<C>>>