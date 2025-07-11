def diagnose_skin_condition():
    """
    Analyzes a clinical case to suggest the most likely diagnosis using a simple scoring model.
    """
    # 1. Define patient's key clinical features from the case study
    patient_profile = {
        "age": 64,
        "BMI": 39,
        "is_smoker": True,
        "lesion_locations": ["axillary folds", "inframammary folds", "inguinal regions"],
        "lesion_types": ["large bullae", "erythematous plaques", "purulent nodules"]
    }

    # 2. Initialize a scoring dictionary for the potential diagnoses
    scores = {
        "Malignant Intertrigo": 0,
        "Allergic contact dermatitis": 0,
        "Hidradenitis Supportiva": 0,
        "Atopic dermatitis": 0,
        "Psoriasis": 0
    }

    # 3. Apply a scoring logic based on clinical knowledge
    # Score Hidradenitis Supportiva (HS)
    # Strong association with obesity and smoking
    if patient_profile["BMI"] > 30:
        scores["Hidradenitis Supportiva"] += 2
    if patient_profile["is_smoker"]:
        scores["Hidradenitis Supportiva"] += 2
    # Classic locations (apocrine gland areas)
    if "axillary folds" in patient_profile["lesion_locations"] and "inguinal regions" in patient_profile["lesion_locations"]:
        scores["Hidradenitis Supportiva"] += 3
    # Hallmark lesion type
    if "purulent nodules" in patient_profile["lesion_types"]:
        scores["Hidradenitis Supportiva"] += 3

    # Score Psoriasis
    # Can present as plaques in intertriginous areas (Inverse Psoriasis)
    if "inframammary folds" in patient_profile["lesion_locations"] and "erythematous plaques" in patient_profile["lesion_types"]:
        scores["Psoriasis"] += 2
    # Purulent nodules are not typical for this type of psoriasis
    if "purulent nodules" in patient_profile["lesion_types"]:
        scores["Psoriasis"] -= 1
        
    # Score Allergic contact dermatitis
    # Bullae can be a feature
    if "large bullae" in patient_profile["lesion_types"]:
        scores["Allergic contact dermatitis"] += 1
    # Purulent nodules are not typical
    if "purulent nodules" in patient_profile["lesion_types"]:
        scores["Allergic contact dermatitis"] -= 2

    # 4. Print the reasoning and the scores ("the final equation")
    print("Clinical Case Analysis:")
    print("=" * 30)
    print(f"Patient BMI: {patient_profile['BMI']} (Obesity is a risk factor for HS)")
    print(f"Patient is a smoker: {patient_profile['is_smoker']} (Smoking is a risk factor for HS)")
    print(f"Lesion Locations: {', '.join(patient_profile['lesion_locations'])} (Classic for HS)")
    print(f"Lesion Types: {', '.join(patient_profile['lesion_types'])} (Purulent nodules are a hallmark of HS)")
    print("=" * 30)
    
    print("\nDiagnostic Score Calculation:")
    print("-" * 30)
    for diagnosis, score in scores.items():
        print(f"Score for {diagnosis}: {score}")
    print("-" * 30)

    # 5. Determine and print the final conclusion
    most_likely_diagnosis = max(scores, key=scores.get)
    print(f"\nConclusion: Based on the scoring, the most likely diagnosis is '{most_likely_diagnosis}'.")

# Execute the function
diagnose_skin_condition()