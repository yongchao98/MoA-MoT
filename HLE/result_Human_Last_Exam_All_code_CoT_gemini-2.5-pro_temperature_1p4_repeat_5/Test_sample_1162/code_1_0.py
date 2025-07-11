import operator

def solve_clinical_case():
    """
    Analyzes a clinical vignette using a scoring system to determine the most likely diagnosis.
    """
    # Patient's key findings
    findings = {
        "Pelvic mass": {"points": 2, "comment": "Common in several pediatric cancers."},
        "Hypertension": {"points": 2, "comment": "Can be caused by renin (Nephroblastoma) or catecholamines (Neuroblastoma)."},
        "Aniridia": {"points": 5, "comment": "A very strong clue for WAGR syndrome, associated with Wilms tumor (Nephroblastoma)."},
        "Developmental delay (speech)": {"points": 1, "comment": "Also a component of WAGR syndrome."},
        "Age (2 years old)": {"points": 0, "comment": "Peak age for both Neuroblastoma and Nephroblastoma."},
        "Exclusion of brain tumor": {"points": -10, "comment": "Pelvic mass makes a primary brain tumor like Astrocytoma extremely unlikely."}
    }

    # Potential diagnoses from the answer choices
    diagnoses = {
        "A. Germ cell tumor": 0,
        "B. Astrocytoma": 0,
        "C. Neuroblastoma": 0,
        "D. Nephroblastoma": 0,
        "E. Ewing sarcoma": 0
    }

    print("Analyzing patient findings with a diagnostic scoring model:\n")
    
    # List to store the numbers for the final equation of the winning diagnosis
    winning_diagnosis_equation_parts = []

    # Score based on Pelvic mass
    mass_points = findings["Pelvic mass"]["points"]
    diagnoses["A. Germ cell tumor"] += mass_points
    diagnoses["C. Neuroblastoma"] += mass_points
    diagnoses["D. Nephroblastoma"] += mass_points
    diagnoses["E. Ewing sarcoma"] += mass_points
    print(f"- Finding: Pelvic mass. Adding {mass_points} points to Germ cell tumor, Neuroblastoma, Nephroblastoma, and Ewing sarcoma.")

    # Score based on Hypertension
    hypertension_points = findings["Hypertension"]["points"]
    diagnoses["C. Neuroblastoma"] += hypertension_points
    diagnoses["D. Nephroblastoma"] += hypertension_points
    print(f"- Finding: Hypertension. Adding {hypertension_points} points to Neuroblastoma and Nephroblastoma.")

    # Score based on Aniridia
    aniridia_points = findings["Aniridia"]["points"]
    diagnoses["D. Nephroblastoma"] += aniridia_points
    print(f"- Finding: Aniridia. This strongly suggests WAGR syndrome. Adding {aniridia_points} points to Nephroblastoma.")

    # Score based on Developmental delay
    delay_points = findings["Developmental delay (speech)"]["points"]
    diagnoses["D. Nephroblastoma"] += delay_points
    print(f"- Finding: Developmental delay. Also associated with WAGR syndrome. Adding {delay_points} points to Nephroblastoma.")
    
    # Score based on exclusion of brain tumor
    exclusion_points = findings["Exclusion of brain tumor"]["points"]
    diagnoses["B. Astrocytoma"] += exclusion_points
    print(f"- Logic: The presence of a pelvic mass rules out a primary brain tumor. Applying {exclusion_points} points to Astrocytoma.")


    # Determine the most likely diagnosis
    most_likely_diagnosis = max(diagnoses.items(), key=operator.itemgetter(1))[0]

    print("\n--- Final Scores ---")
    for diagnosis, score in diagnoses.items():
        print(f"{diagnosis}: {score} points")
    print("--------------------\n")

    # Construct and print the final equation for the winning diagnosis
    print("Conclusion: The highest score indicates the most likely diagnosis.")
    
    # Build the equation string for the winner
    total_score_for_winner = diagnoses[most_likely_diagnosis]
    equation_str = f"Final equation for Nephroblastoma: {mass_points} (from mass) + {hypertension_points} (from hypertension) + {aniridia_points} (from aniridia) + {delay_points} (from delay) = {total_score_for_winner}"
    print(equation_str)
    
    print(f"\nThe patient's presentation of aniridia, developmental delay, hypertension, and a pelvic mass is classic for WAGR syndrome, which includes Wilms tumor (Nephroblastoma).")
    print(f"The most likely diagnosis is {most_likely_diagnosis}.")


solve_clinical_case()