import sys

def evaluate_diagnoses():
    """
    Analyzes clinical findings to score potential diagnoses.
    This function simulates a diagnostic process by assigning points to each diagnosis
    based on how well it matches the key findings from the case study.
    """

    # --- Key Clinical Findings from the Case ---
    # We assign weights based on the diagnostic significance of each finding.
    findings = {
        "Difficult Colonoscopy": {"points": 1, "present": True},
        "LUQ & Left Shoulder Pain": {"points": 5, "present": True}, # Highly specific (Kehr's Sign)
        "Hemorrhagic Shock": {"points": 3, "present": True}, # Indicates severe, acute bleeding
        "Peritoneal Signs": {"points": 2, "present": True},
        "No Polypectomy Performed": {"points": -100, "present": True} # Exclusionary finding
    }

    # --- Potential Diagnoses ---
    diagnoses = {
        "A. Colonic perforation": {
            "explanation": "Pain location is atypical; massive bleeding less common.",
            "score_components": ["Difficult Colonoscopy", "Peritoneal Signs"]
        },
        "B. Lower GI bleeding": {
            "explanation": "Pain location points to an intraperitoneal, not intraluminal, source.",
            "score_components": ["Hemorrhagic Shock"],
            "negative_components": ["LUQ & Left Shoulder Pain"]
        },
        "C. Splenic laceration": {
            "explanation": "Classic presentation: difficult colonoscopy, LUQ/shoulder pain, and shock.",
            "score_components": ["Difficult Colonoscopy", "LUQ & Left Shoulder Pain", "Hemorrhagic Shock", "Peritoneal Signs"]
        },
        "D. Postpolypectomy syndrome": {
            "explanation": "Ruled out by the absence of a polypectomy.",
            "score_components": ["No Polypectomy Performed"]
        }
    }

    print("Evaluating diagnoses based on clinical findings...\n")
    
    final_scores = {}

    for diagnosis, details in diagnoses.items():
        total_score = 0
        equation_str = f"{diagnosis}: "
        
        # Add points for matching findings
        for finding_name in details.get("score_components", []):
            if findings[finding_name]["present"]:
                points = findings[finding_name]["points"]
                total_score += points
                equation_str += f"{points} (for {finding_name}) + "

        # Subtract points for contradictory findings
        for finding_name in details.get("negative_components", []):
             if findings[finding_name]["present"]:
                # The pain location argues against this diagnosis, so we subtract points
                points = findings[finding_name]["points"]
                total_score -= points
                equation_str += f" -{points} (against due to {finding_name}) + "


        # Clean up the equation string and print
        if equation_str.endswith(" + "):
            equation_str = equation_str[:-3]
        equation_str += f" = {total_score}"
        
        print(equation_str)
        print(f"   Rationale: {details['explanation']}\n")
        final_scores[diagnosis] = total_score

    # Determine the most likely diagnosis
    most_likely_diagnosis = max(final_scores, key=final_scores.get)
    
    print("-" * 40)
    print(f"The most likely diagnosis is {most_likely_diagnosis} with a score of {final_scores[most_likely_diagnosis]}.")
    print("-" * 40)


if __name__ == "__main__":
    evaluate_diagnoses()
    # The final answer is determined by the script's output.
    # The clinical presentation is a textbook case for iatrogenic splenic laceration
    # following a colonoscopy, making C the correct choice.
    sys.stdout.flush() # Ensure all prints are displayed before the final answer tag.
    # The final answer is represented by the letter corresponding to the highest-scoring diagnosis.
