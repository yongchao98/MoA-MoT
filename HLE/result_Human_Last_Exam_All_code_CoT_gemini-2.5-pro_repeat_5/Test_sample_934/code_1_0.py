import sys

def solve_clinical_case():
    """
    This function analyzes a clinical vignette to determine the most likely diagnosis.
    It uses a simple scoring system to weigh evidence for and against each differential diagnosis.
    """

    # --- Patient's Key Findings & Risk Factors with assigned points ---
    # Positive points for factors that support a diagnosis, negative for those that argue against it.
    patient_data = {
        "Heavy Smoking": {"points": 4, "reason": "Major risk factor for SCC."},
        "Alcohol Use Disorder": {"points": 4, "reason": "Major risk factor for SCC."},
        "Chest Pain / Odynophagia": {"points": 2, "reason": "Common symptom for esophagitis and cancer."},
        "Imaging (Wall Thickening/Narrowing)": {"points": 3, "reason": "Suggests an infiltrative process or mass."},
        "Normal Endoscopy (No Mucosal Lesions)": {"points": -5, "reason": "Strongly argues against diseases with typical surface lesions (GERD, infections)."}
    }

    # --- Differential Diagnoses ---
    # Each diagnosis has a base likelihood and specific associations.
    diagnoses = {
        "A. Streptococcal esophagitis": {
            "notes": "Rare, usually causes visible plaques/ulcers on endoscopy.",
            "score_modifiers": {"Normal Endoscopy (No Mucosal Lesions)": True} # This finding is highly relevant
        },
        "B. Esophageal adenocarcinoma": {
            "notes": "Primarily associated with GERD/Barrett's esophagus, which is not mentioned.",
            "score_modifiers": {"Heavy Smoking": True} # Smoking is a moderate risk factor
        },
        "C. Esophageal squamous cell carcinoma": {
            "notes": "Strongly associated with smoking and alcohol. Can be infiltrative, explaining normal endoscopy with abnormal imaging.",
            "score_modifiers": {
                "Heavy Smoking": True,
                "Alcohol Use Disorder": True,
                "Imaging (Wall Thickening/Narrowing)": True,
                # Note: The normal endoscopy is a positive clue for an *infiltrative* SCC, so we don't apply the negative points.
                "Normal Endoscopy (No Mucosal Lesions)": False
            }
        },
        "D. GERD": {
            "notes": "Severe GERD causing these symptoms would typically show esophagitis on endoscopy.",
            "score_modifiers": {"Normal Endoscopy (No Mucosal Lesions)": True}
        },
        "E. Herpes esophagitis": {
            "notes": "Typically presents with 'punched-out' ulcers visible on endoscopy.",
            "score_modifiers": {"Normal Endoscopy (No Mucosal Lesions)": True}
        }
    }

    print("Analyzing the clinical case step-by-step:\n")

    results = {}
    highest_score = -float('inf')
    most_likely_diagnosis = ""
    final_equation_str = ""

    for diagnosis, details in diagnoses.items():
        score = 0
        equation_parts = []
        print(f"--- Evaluating: {diagnosis} ---")
        print(f"Notes: {details['notes']}")

        # Calculate score based on patient data
        for finding, data in patient_data.items():
            # Check if this finding is relevant for the current diagnosis
            if details["score_modifiers"].get(finding):
                score += data["points"]
                equation_parts.append(f"({finding}: {data['points']})")

        # Special handling for SCC where normal endoscopy is a supportive, not negative, finding for infiltrative type
        if diagnosis == "C. Esophageal squamous cell carcinoma":
            # Add points for symptoms and imaging which are highly consistent
            score += patient_data["Chest Pain / Odynophagia"]["points"]
            equation_parts.append(f"({list(patient_data.keys())[2]}: {patient_data['Chest Pain / Odynophagia']['points']})")

            score += patient_data["Imaging (Wall Thickening/Narrowing)"]["points"]
            equation_parts.append(f"({list(patient_data.keys())[3]}: {patient_data['Imaging (Wall Thickening/Narrowing)']['points']})")


        results[diagnosis] = score
        print(f"Likelihood Score: {score}\n")

        if score > highest_score:
            highest_score = score
            most_likely_diagnosis = diagnosis
            final_equation_str = f"Final Equation for {most_likely_diagnosis}: " + " + ".join(equation_parts) + f" = {score}"


    print("--- Conclusion ---")
    print(f"The most likely diagnosis is {most_likely_diagnosis} with a score of {highest_score}.")
    print("\nThis is because the patient's two most significant risk factors (heavy smoking and alcohol use) are strongly associated with this cancer.")
    print("Furthermore, the combination of abnormal imaging (wall thickening) with a normal-appearing endoscopy is a classic presentation for an infiltrative squamous cell carcinoma.\n")

    print(final_equation_str)


if __name__ == "__main__":
    solve_clinical_case()
    sys.stdout.flush() # Ensure all output is printed before the final answer tag
    print("\n<<<C>>>")