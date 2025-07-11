def solve_diagnosis():
    """
    Analyzes clinical findings to determine the most likely diagnosis
    by scoring each option based on key evidence from the case.
    """
    # Key clinical findings from the case text.
    # A value of 1 represents presence, -1 represents a direct contradiction.
    findings = {
        'No Polypectomy': -1,  # A crucial negative finding
        'LUQ Pain': 1,
        'Left Shoulder Pain (Kehr\'s Sign)': 1,
        'Rapid Hemoglobin Drop & Shock': 1,
        'Peritoneal Signs': 1
    }

    # Diagnoses and their evidence mapping.
    # 1: a key supporting feature
    # -1: a key contradicting feature
    # 0: neutral or less likely
    diagnoses_evidence_map = {
        'A. Colonic perforation': {
            'No Polypectomy': 0,
            'LUQ Pain': 1,  # Possible if at splenic flexure
            'Left Shoulder Pain (Kehr\'s Sign)': -1,  # More classic for blood than air/stool
            'Rapid Hemoglobin Drop & Shock': -1, # Bleeding of this magnitude is less typical
            'Peritoneal Signs': 1
        },
        'B. Lower GI bleeding': {
            'No Polypectomy': -1, # No source of bleed mentioned
            'LUQ Pain': -1, # Peritoneal signs are not typical for intraluminal bleeding
            'Left Shoulder Pain (Kehr\'s Sign)': -1,
            'Rapid Hemoglobin Drop & Shock': 1,
            'Peritoneal Signs': -1
        },
        'C. Splenic laceration': {
            'No Polypectomy': 0, # Unrelated to mechanism
            'LUQ Pain': 1, # Classic location
            'Left Shoulder Pain (Kehr\'s Sign)': 1, # Classic sign of diaphragmatic irritation from blood
            'Rapid Hemoglobin Drop & Shock': 1, # Classic for massive hemorrhage from spleen
            'Peritoneal Signs': 1 # Classic for hemoperitoneum
        },
        'D. Postpolypectomy syndrome': {
            'No Polypectomy': -1, # This finding EXCLUDES the diagnosis
            'LUQ Pain': 0,
            'Left Shoulder Pain (Kehr\'s Sign)': -1,
            'Rapid Hemoglobin Drop & Shock': -1, # Not a hemorrhagic syndrome
            'Peritoneal Signs': 1 # Inflammation can cause this, but it's ruled out
        }
    }

    scores = {}
    print("Calculating Likelihood Score for Each Diagnosis:\n")

    for diagnosis, evidence in diagnoses_evidence_map.items():
        score = 0
        equation_parts = []
        for finding_name, finding_present in findings.items():
            # For the special case of 'No Polypectomy', we just use its evidence score directly.
            if finding_name == 'No Polypectomy':
                weight = evidence[finding_name]
            else:
                # For other findings, the weight is the finding's presence (1) * its relevance (evidence score).
                weight = finding_present * evidence[finding_name]

            score += weight
            # Append the number for the equation string
            equation_parts.append(str(weight))

        scores[diagnosis] = score
        # The print statement fulfills the "output each number in the final equation" requirement
        print(f"Diagnosis: {diagnosis}")
        print(f"  Score Calculation: {' + '.join(equation_parts).replace('+ -', '- ')} = {score}")
        print("-" * 20)


    # Find the diagnosis with the highest score
    most_likely_diagnosis = max(scores, key=scores.get)
    
    # Extract the letter from the diagnosis string
    final_answer_letter = most_likely_diagnosis.split('.')[0]

    print(f"\nThe most likely diagnosis based on the scoring is: {most_likely_diagnosis}")
    print(f"<<<{final_answer_letter}>>>")

solve_diagnosis()