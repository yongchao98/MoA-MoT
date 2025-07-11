def solve_diagnosis():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring potential diagnoses against patient findings.
    """
    # Step 1: Define patient findings from the clinical case.
    # BMI 39 is obesity. The patient is a smoker. Locations and lesions are noted.
    patient_findings = {
        "locations": {"axillary", "inframammary", "inguinal"},
        "lesions": {"bullae", "erythematous plaques", "purulent nodules"},
        "risk_factors": {"obesity", "smoking"}
    }

    # Step 2: Define diagnostic criteria for each condition.
    # Key features are highly suggestive, supporting features are common but less specific,
    # and contradictions argue against the diagnosis.
    diagnoses_criteria = {
        "A. Malignant Intertrigo": {
            "key_lesions": {"persistent_plaques"},
            "key_locations": {"inframammary", "inguinal", "axillary"},
            "supporting_risk_factors": set(),
            "contradictions": {"purulent nodules"} # Not a typical feature
        },
        "B. Allergic contact dermatitis": {
            "key_lesions": {"vesicles", "eczema"},
            "key_locations": set(), # Variable location
            "supporting_risk_factors": set(),
            "contradictions": {"purulent nodules"} # Not a feature
        },
        "C. Hidradenitis Supportiva": {
            "key_lesions": {"purulent nodules"},
            "key_locations": {"axillary", "inguinal", "inframammary"},
            "supporting_risk_factors": {"obesity", "smoking"}
        },
        "D. Atopic dermatitis": {
            "key_lesions": {"eczemateous_rash"},
            "key_locations": {"flexural"},
            "supporting_risk_factors": set(),
            "contradictions": {"purulent nodules"} # Not a feature
        },
        "E. Psoriasis": {
            "key_lesions": {"erythematous plaques"}, # Inverse Psoriasis
            "key_locations": {"axillary", "inguinal", "inframammary"},
            "supporting_risk_factors": set(),
            "contradictions": {"purulent nodules"} # Not a feature of Inverse Psoriasis
        }
    }

    # Step 3: Calculate a score for each diagnosis and print the reasoning.
    print("Analyzing clinical case based on key features...")
    print("Scoring model: (Matching Key Lesions * 3) + (Matching Key Locations * 2) + (Matching Risk Factors * 1) - (Contradictions * 5)")
    print("-" * 110)

    final_scores = {}
    for diagnosis, rules in diagnoses_criteria.items():
        # Score matching key lesions
        lesion_matches = patient_findings["lesions"].intersection(rules.get("key_lesions", set()))
        lesion_score = len(lesion_matches) * 3

        # Score matching key locations
        location_matches = patient_findings["locations"].intersection(rules.get("key_locations", set()))
        location_score = len(location_matches) * 2

        # Score matching risk factors
        risk_matches = patient_findings["risk_factors"].intersection(rules.get("supporting_risk_factors", set()))
        risk_score = len(risk_matches) * 1
        
        # Apply penalty for contradictory findings
        contradiction_matches = patient_findings["lesions"].intersection(rules.get("contradictions", set()))
        penalty = len(contradiction_matches) * 5

        total_score = lesion_score + location_score + risk_score - penalty
        final_scores[diagnosis] = total_score
        
        # Print the "equation" for each diagnosis, fulfilling the user's request
        print(f"Diagnosis: {diagnosis:<30} | Score Equation: {lesion_score} + {location_score} + {risk_score} - {penalty} = {total_score}")

    # Step 4: Determine the best diagnosis.
    best_diagnosis = max(final_scores, key=final_scores.get)

    print("-" * 110)
    print("\nConclusion:")
    print("The patient's presentation of purulent nodules in classic intertriginous sites (axilla, inguinal, inframammary)")
    print("combined with major risk factors (obesity and smoking) strongly points to a specific diagnosis.")
    print(f"\nThe diagnosis with the highest score is: {best_diagnosis}")

solve_diagnosis()