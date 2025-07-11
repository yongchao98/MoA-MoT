def diagnose_skin_condition():
    """
    Analyzes a patient's clinical data to suggest the most likely dermatological diagnosis.
    This function uses a scoring system to evaluate how patient findings match the criteria
    for several potential diagnoses.
    """

    # 1. Define Patient's Profile from the case study
    patient_profile = {
        "lesions": ["large bullae", "erythematous plaques", "purulent nodules"],
        "locations": ["axillary folds", "inframammary folds", "inguinal regions"],
        "risk_factors": ["obesity", "smoking", "type 2 diabetes"]
    }

    # 2. Define Diagnostic Criteria with scoring points.
    # Key features receive higher scores. Negative features receive penalties.
    diagnostic_criteria = {
        "Malignant Intertrigo": {"scores": {"plaques": 1, "history_of_cancer": 2}, "mismatch_penalty": {"purulent nodules": -3}},
        "Allergic contact dermatitis": {"scores": {"plaques": 1, "bullae": 1}, "mismatch_penalty": {"purulent nodules": -3}},
        "Hidradenitis Supportiva": {"scores": {"purulent nodules": 3, "intertriginous_locations": 3, "obesity": 2, "smoking": 2, "plaques": 1}, "mismatch_penalty": {}},
        "Atopic dermatitis": {"scores": {"plaques": 1}, "mismatch_penalty": {"purulent nodules": -2, "bullae": -2}},
        "Psoriasis": {"scores": {"plaques": 2, "intertriginous_locations": 1}, "mismatch_penalty": {"purulent nodules": -3, "bullae": -2}},
    }
    
    print("Analyzing patient data to determine the most likely diagnosis...\n")
    final_scores = {}

    # 3. Calculate score for each diagnosis
    for diagnosis, criteria in diagnostic_criteria.items():
        score = 0
        score_details = []
        
        # Match lesions
        if "purulent nodules" in patient_profile["lesions"] and "purulent nodules" in criteria["scores"]:
            points = criteria["scores"]["purulent nodules"]
            score += points
            score_details.append(str(points))
        if "large bullae" in patient_profile["lesions"] and "bullae" in criteria["scores"]:
            points = criteria["scores"]["bullae"]
            score += points
            score_details.append(str(points))
        if "erythematous plaques" in patient_profile["lesions"] and "plaques" in criteria["scores"]:
            points = criteria["scores"]["plaques"]
            score += points
            score_details.append(str(points))

        # Match locations
        if "intertriginous_locations" in criteria["scores"]:
            points = criteria["scores"]["intertriginous_locations"]
            score += points
            score_details.append(str(points))
            
        # Match risk factors
        if "obesity" in patient_profile["risk_factors"] and "obesity" in criteria["scores"]:
            points = criteria["scores"]["obesity"]
            score += points
            score_details.append(str(points))
        if "smoking" in patient_profile["risk_factors"] and "smoking" in criteria["scores"]:
            points = criteria["scores"]["smoking"]
            score += points
            score_details.append(str(points))
            
        # Apply penalties for mismatched features
        if "purulent nodules" in patient_profile["lesions"] and "purulent nodules" in criteria.get("mismatch_penalty", {}):
            penalty = criteria["mismatch_penalty"]["purulent nodules"]
            score += penalty
            score_details.append(str(penalty))
        if "large bullae" in patient_profile["lesions"] and "bullae" in criteria.get("mismatch_penalty", {}):
            penalty = criteria["mismatch_penalty"]["bullae"]
            score += penalty
            score_details.append(str(penalty))
            
        final_scores[diagnosis] = score
        
        # 4. Print the reasoning and final score for each diagnosis
        print(f"--- {diagnosis} ---")
        if not score_details:
             print("No matching features found.")
             print("Final Score: 0")
        else:
            equation = " + ".join(score_details).replace("+ -", "- ")
            print(f"Scoring Equation: {equation} = {score}")
            print(f"Final Score: {score}")
        print("-" * (len(diagnosis) + 6) + "\n")

    # 5. Determine and print the most likely diagnosis
    best_diagnosis = max(final_scores, key=final_scores.get)
    print(f"\nConclusion: The diagnosis with the highest score is '{best_diagnosis}'.")

diagnose_skin_condition()