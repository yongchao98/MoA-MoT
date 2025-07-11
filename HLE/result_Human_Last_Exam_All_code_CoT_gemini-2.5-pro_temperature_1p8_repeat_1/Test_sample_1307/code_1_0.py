import sys

def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by matching patient symptoms to known medical conditions.
    """
    # Step 1: Define patient's key findings from the vignette
    patient_findings = {
        "palpable thrill (vibration)",
        "continuous bruit (murmur)",
        "recent femoral artery access"
    }

    # Step 2: Define the differential diagnoses and their characteristic signs
    diagnoses = {
        'A': {
            "name": "Femoral venous thrombosis",
            "symptoms": {"leg swelling (edema)", "pain", "warmth"}
        },
        'B': {
            "name": "Arterial embolism",
            "symptoms": {"acute limb ischemia", "pulselessness", "pallor"}
        },
        'C': {
            "name": "Retroperitoneal hematoma",
            "symptoms": {"flank pain", "hypotension", "anemia"}
        },
        'D': {
            "name": "Femoral artery dissection",
            "symptoms": {"acute limb ischemia", "severe pain"}
        },
        'E': {
            "name": "Hamartoma",
            "symptoms": {"developmental benign mass"}
        },
        'F': {
            "name": "Femoral artery pseudoaneurysm",
            "symptoms": {"palpable thrill (vibration)", "continuous bruit (murmur)", "pulsatile mass", "recent femoral artery access"}
        },
        'G': {
            "name": "None of these choices",
            "symptoms": set()
        },
        'H': {
            "name": "Arterio-capillary communication",
            "symptoms": {"microvascular changes"}
        }
    }

    # Step 3 & 4: Score each diagnosis based on matching findings
    best_match_choice = None
    max_score = -1
    best_match_details = {}

    for choice, details in diagnoses.items():
        # Find the intersection of patient findings and diagnosis symptoms
        matching_symptoms = patient_findings.intersection(details["symptoms"])
        score = len(matching_symptoms)

        if score > max_score:
            max_score = score
            best_match_choice = choice
            best_match_details = {
                "name": details["name"],
                "matching_symptoms": list(matching_symptoms),
                "score": score
            }

    # Step 5: Output the reasoning and the result
    print("Analyzing patient findings against possible diagnoses...")
    print(f"Patient's key findings: {', '.join(patient_findings)}")
    print("-" * 30)

    print(f"Best match found: {best_match_details['name']}")
    print("This diagnosis is selected because it accounts for the following patient findings:")
    for symptom in best_match_details['matching_symptoms']:
        print(f"- {symptom}")

    print("\nFinal score calculation for the best match:")
    equation_parts = []
    for i in range(best_match_details['score']):
        equation_parts.append("1")
    
    if not equation_parts:
        print("Score = 0")
    else:
        # This loop prints each number in the final equation
        for i, part in enumerate(equation_parts):
            sys.stdout.write(part)
            if i < len(equation_parts) - 1:
                sys.stdout.write(" + ")
        print(f" = {best_match_details['score']}")

    print("\nConclusion: The combination of a palpable thrill and a continuous bruit after arterial puncture is the classic presentation of a pseudoaneurysm.")
    print(f"The correct choice is {best_match_choice}.")


solve_clinical_case()