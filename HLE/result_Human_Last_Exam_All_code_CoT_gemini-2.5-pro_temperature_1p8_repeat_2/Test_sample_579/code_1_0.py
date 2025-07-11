import collections

def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """

    # 1. Define patient's key features from the vignette
    patient_profile = {
        "locations": ["axillary", "inframammary", "inguinal"],
        "lesions": ["bullae/abscesses", "plaques", "purulent nodules"],
        "risk_factors": ["obesity", "smoking"]
    }

    # 2. Define profiles for each potential diagnosis
    diagnostic_profiles = {
        "Malignant Intertrigo": {
            "locations": ["inframammary", "inguinal", "axillary"],
            "lesions": ["plaques", "nodules", "ulceration"],
            "risk_factors": ["history of cancer"],
            "notes": "Less likely to cause purulent nodules or involve all three sites simultaneously as a primary diagnosis."
        },
        "Allergic contact dermatitis": {
            "locations": ["any"],
            "lesions": ["erythema", "vesicles", "plaques"],
            "risk_factors": ["allergen exposure"],
            "notes": "Purulent nodules and deep abscesses are not typical."
        },
        "Hidradenitis Suppurativa": {
            "locations": ["axillary", "inframammary", "inguinal"],
            "lesions": ["purulent nodules", "abscesses", "sinus tracts", "scarring"],
            "risk_factors": ["obesity", "smoking", "female sex"],
            "notes": "Classic presentation involves intertriginous sites with purulent, painful nodules and abscesses."
        },
        "Atopic dermatitis": {
            "locations": ["flexural (e.g., antecubital)"],
            "lesions": ["plaques", "lichenification", "pruritus"],
            "risk_factors": ["personal/family history of atopy"],
            "notes": "Purulent nodules and bullae are not primary features."
        },
        "Psoriasis (Inverse)": {
            "locations": ["axillary", "inframammary", "inguinal"],
            "lesions": ["plaques (well-demarcated, erythematous, non-scaly)"],
            "risk_factors": ["family history", "autoimmune"],
            "notes": "Does not typically present with purulent nodules or bullae/abscesses."
        }
    }

    scores = collections.defaultdict(int)
    analysis_log = []

    print("Analyzing patient profile against potential diagnoses...\n")

    # 3 & 4. Score each diagnosis based on matching features
    for diagnosis, profile in diagnostic_profiles.items():
        analysis_log.append(f"--- Evaluating: {diagnosis} ---")
        
        # Score locations
        location_matches = set(patient_profile["locations"]) & set(profile["locations"])
        if profile["locations"] == ["any"]: # For non-specific locations like ACD
             location_matches.add("site match")
        score = len(location_matches)
        analysis_log.append(f"Location matches: {len(location_matches)} ({', '.join(location_matches) if location_matches else 'None'})")

        # Score lesions
        lesion_matches = set(patient_profile["lesions"]) & set(profile["lesions"])
        # HS abscesses can be described as bullae
        if diagnosis == "Hidradenitis Suppurativa" and "bullae/abscesses" in patient_profile["lesions"]:
            lesion_matches.add("bullae/abscesses")
        score += len(lesion_matches)
        analysis_log.append(f"Lesion matches: {len(lesion_matches)} ({', '.join(lesion_matches) if lesion_matches else 'None'})")
        
        # Score risk factors
        risk_matches = set(patient_profile["risk_factors"]) & set(profile["risk_factors"])
        score += len(risk_matches)
        analysis_log.append(f"Risk factor matches: {len(risk_matches)} ({', '.join(risk_matches) if risk_matches else 'None'})")
        
        scores[diagnosis] = score
        analysis_log.append(f"Note: {profile['notes']}")
        analysis_log.append(f"Total Score for {diagnosis}: {score}\n")

    # 5. Print the step-by-step analysis
    for line in analysis_log:
        print(line)

    # 6. Determine and print the final answer
    most_likely_diagnosis = max(scores, key=scores.get)

    print("="*30)
    print("Conclusion:")
    print("The patient presents with lesions (purulent nodules, plaques, large bullae/abscesses) in classic intertriginous areas (axillary, inframammary, inguinal).")
    print("The patient also has major risk factors: obesity (BMI 39) and smoking.")
    print("This combination of findings is the hallmark of Hidradenitis Suppurativa.")
    print(f"\nThe diagnosis with the highest match score is: {most_likely_diagnosis}")
    print("="*30)

solve_medical_case()
<<<C>>>