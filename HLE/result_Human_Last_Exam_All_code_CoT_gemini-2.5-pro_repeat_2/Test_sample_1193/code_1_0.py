def evaluate_diagnoses():
    """
    Evaluates potential diagnoses for a post-operative patient with hypoxemia
    by scoring them against key clinical findings.
    """

    # Clinical findings from the case vignette
    # A score of 10 indicates a strong match, 0 indicates a poor match or contradiction.
    timeline_is_delayed_post_op = {"score": 10, "reason": "Timeline (29 days)"}
    is_high_risk_surgery = {"score": 10, "reason": "High-Risk Surgery (Whipple)"}
    symptoms_match_ards = {"score": 10, "reason": "Symptoms (Hypoxemia + Crackles)"}
    timeline_is_acute = {"score": 0, "reason": "Timeline (29 days)"} # Contradicts acute events

    diagnoses = {
        "A. Acute blood transfusion reaction": [timeline_is_acute, is_high_risk_surgery, symptoms_match_ards],
        "B. Iodine-related reaction": [timeline_is_acute, is_high_risk_surgery, symptoms_match_ards],
        "C. Sensitivity reaction": [timeline_is_delayed_post_op, is_high_risk_surgery, symptoms_match_ards],
        "D. Sepsis": [timeline_is_delayed_post_op, is_high_risk_surgery, symptoms_match_ards],
        "E. Myocyte necrosis": [timeline_is_delayed_post_op, is_high_risk_surgery, symptoms_match_ards],
        "F. Respiratory deconditioning": [timeline_is_delayed_post_op, is_high_risk_surgery, symptoms_match_ards],
        "G. Lung exhaustion": [timeline_is_delayed_post_op, is_high_risk_surgery, symptoms_match_ards],
        "H. Air pollution sensitivity": [timeline_is_delayed_post_op, is_high_risk_surgery, symptoms_match_ards],
    }

    # Assigning specific weights for each diagnosis based on medical knowledge
    weights = {
        "A. Acute blood transfusion reaction": [0, 0, 1], # Timeline is wrong
        "B. Iodine-related reaction": [0, 0, 1], # Timeline is wrong
        "C. Sensitivity reaction": [0.5, 0.2, 0.5], # Vague, less likely to cause ARDS
        "D. Sepsis": [1, 1, 1], # Perfect match for all findings
        "E. Myocyte necrosis": [0.2, 0.2, 0.3], # Doesn't typically present this way
        "F. Respiratory deconditioning": [1, 0.5, 0.1], # Doesn't cause acute ARDS
        "G. Lung exhaustion": [1, 0.5, 0.1], # Not a real diagnosis, similar to F
        "H. Air pollution sensitivity": [0.1, 0, 0.2], # Unlikely to cause this in post-op setting
    }

    best_diagnosis = ""
    max_score = -1

    print("Evaluating diagnoses based on clinical evidence...\n")

    for diagnosis, factors in diagnoses.items():
        score = 0
        equation_str = f"Diagnosis: {diagnosis}\nScore = "
        
        # Unpack the three main factors for scoring
        timeline_factor = factors[0]
        surgery_factor = factors[1]
        symptoms_factor = factors[2]

        # Get the specific weights for this diagnosis
        w = weights[diagnosis]

        # Calculate scores for each component
        timeline_score = timeline_factor["score"] * w[0]
        surgery_score = surgery_factor["score"] * w[1]
        symptoms_score = symptoms_factor["score"] * w[2]
        
        # Final score is the sum of weighted component scores
        total_score = timeline_score + surgery_score + symptoms_score

        # Build the equation string showing the numbers
        equation_str += f"({timeline_factor['reason']}: {timeline_score:.1f}) + "
        equation_str += f"({surgery_factor['reason']}: {surgery_score:.1f}) + "
        equation_str += f"({symptoms_factor['reason']}: {symptoms_score:.1f})"
        
        print(equation_str)
        print(f"Final Score = {total_score:.1f}\n" + "-"*40)

        if total_score > max_score:
            max_score = total_score
            best_diagnosis = diagnosis

    print(f"\nConclusion: The most likely diagnosis is '{best_diagnosis}' with a score of {max_score:.1f}.")

evaluate_diagnoses()