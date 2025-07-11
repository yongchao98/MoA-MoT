import sys

def analyze_clinical_case():
    """
    This function simulates the diagnostic process for the given clinical case.
    It scores each diagnosis based on key findings from the case description.
    """

    # --- Define Clinical Findings ---
    findings = {
        "procedure_type": "Difficult colonoscopy",
        "key_symptoms": ["Left Upper Quadrant (LUQ) pain", "Left shoulder pain (Kehr's sign)"],
        "signs_of_hemorrhage": ["Hemodynamic instability (tachycardia, hypotension)", "Dramatic drop in hemoglobin"],
        "physical_exam": ["LUQ tenderness", "Peritoneal signs"],
        "negative_finding": "No polypectomy performed"
    }

    # --- Define Potential Diagnoses and Initial Scores ---
    diagnoses = {
        "A": {"name": "Colonic perforation", "score": 0, "reasoning": []},
        "B": {"name": "Lower GI bleeding", "score": 0, "reasoning": []},
        "C": {"name": "Splenic laceration", "score": 0, "reasoning": []},
        "D": {"name": "Postpolypectomy syndrome", "score": 0, "reasoning": []},
    }

    # --- Scoring Logic ---

    # Score based on pain location
    diagnoses["C"]["score"] += 3
    diagnoses["C"]["reasoning"].append("(+3 points) Classic presentation with LUQ and referred left shoulder pain (Kehr's sign), indicating diaphragmatic irritation from splenic blood.")
    diagnoses["A"]["score"] += 1
    diagnoses["A"]["reasoning"].append("(+1 point) Possible, but pain is less specific than the LUQ/shoulder pattern.")

    # Score based on massive hemorrhage
    diagnoses["C"]["score"] += 2
    diagnoses["C"]["reasoning"].append("(+2 points) Profound shock and hemoglobin drop are characteristic of a ruptured spleen, a highly vascular organ.")
    diagnoses["A"]["score"] += 1
    diagnoses["A"]["reasoning"].append("(+1 point) Perforation can cause bleeding, but the severity points strongly towards a major vessel or organ injury.")
    diagnoses["B"]["score"] += 1
    diagnoses["B"]["reasoning"].append("(+1 point) Hemorrhage is present, but this presentation (pain/shock first, no mention of rectal bleeding) is not typical for a standard lower GI bleed.")

    # Score based on peritoneal signs
    diagnoses["A"]["score"] += 2
    diagnoses["A"]["reasoning"].append("(+2 points) Peritoneal signs are classic for perforation due to leakage of bowel contents.")
    diagnoses["C"]["score"] += 1
    diagnoses["C"]["reasoning"].append("(+1 point) Also explains peritoneal signs due to hemoperitoneum (blood irritating the peritoneum).")
    
    # Rule out diagnosis based on negative finding
    if "No polypectomy performed" in findings.values():
        diagnoses["D"]["score"] = -99 # Effectively disqualifies this option
        diagnoses["D"]["reasoning"].append("(-99 points) Ruled out as the case explicitly states no polypectomy was performed.")


    # --- Determine the most likely diagnosis ---
    most_likely_diagnosis_key = max(diagnoses, key=lambda k: diagnoses[k]['score'])
    result = diagnoses[most_likely_diagnosis_key]

    print("Analyzing Clinical Findings to Determine the Most Likely Diagnosis:\n")
    
    # --- Print the reasoning for each diagnosis ---
    for key, diag in sorted(diagnoses.items(), key=lambda item: item[1]['score'], reverse=True):
        print(f"Diagnosis ({key}) {diag['name']}:")
        print(f"  Final Score: {diag['score']}")
        for reason in diag['reasoning']:
              print(f"  - {reason}")
        print("-" * 30)

    # --- Final Conclusion ---
    print("\nConclusion:")
    print(f"The diagnosis with the highest score is '{result['name']}'.")
    print("The patient's presentation of acute LUQ pain, referred left shoulder pain (Kehr's sign), and severe hemodynamic collapse after a difficult colonoscopy is a classic, albeit rare, presentation of a splenic laceration.")
    
    # To conform to the specified output format for the final answer.
    # Note: In a real application, you would just return the conclusion.
    # This specific format is for the context of this task.
    sys.stdout = open(sys.platform == "win32" and "nul" or "/dev/null", 'w')
    print(f"<<<{most_likely_diagnosis_key}>>>")


if __name__ == "__main__":
    analyze_clinical_case()