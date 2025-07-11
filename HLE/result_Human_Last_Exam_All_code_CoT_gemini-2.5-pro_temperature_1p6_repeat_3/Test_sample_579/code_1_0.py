def diagnose_skin_condition():
    """
    Analyzes a clinical case to provide a diagnosis using a scoring system.
    """
    # Step 1: Define key clinical findings from the case.
    patient_info = {
        "BMI": 39,
        "Smokes": True,
        "Locations": ["axillary folds", "inframammary folds", "inguinal regions"],
        "Lesions": ["large bullae", "erythematous plaques", "purulent nodules"]
    }

    # Step 2: Explain the reasoning and scoring system.
    print("Analyzing the clinical case to determine the diagnosis...")
    print("The diagnosis can be determined by evaluating key features. We will assign a point for each major category of findings consistent with Hidradenitis Suppurativa (HS).\n")

    # Score for characteristic locations (axilla, inguinal, inframammary)
    # The patient presents with lesions in all three classic intertriginous areas.
    location_score = 1
    print(f"Finding 1: Lesions in classic HS locations (axilla, groin, inframammary).")
    print(f"Result: Positive. Score = {location_score}\n")

    # Score for characteristic lesions (purulent nodules are a hallmark of HS)
    # The patient has purulent nodules.
    lesion_score = 1
    print(f"Finding 2: Presence of typical HS lesions (e.g., inflammatory/purulent nodules).")
    print(f"Result: Positive. Score = {lesion_score}\n")
    
    # Score for major risk factors (obesity and smoking are strongly associated with HS)
    # The patient has a BMI of 39 and is a smoker.
    risk_factor_score = 1
    print(f"Finding 3: Presence of major HS risk factors (Obesity and Smoking).")
    print(f"Result: Positive. Score = {risk_factor_score}\n")

    # Step 3: Display the "final equation" to show how the evidence combines.
    total_score = location_score + lesion_score + risk_factor_score
    print("The diagnostic evidence can be summarized with the following equation:")
    print(f"{location_score} (Locations) + {lesion_score} (Lesions) + {risk_factor_score} (Risk Factors) = {total_score} (Total Score for HS)")

    # Step 4: Provide the final conclusion based on the evidence.
    print("\nConclusion: The patient's presentation with purulent nodules in the axillary, inguinal, and inframammary folds, combined with risk factors of obesity and smoking, is a classic picture of Hidradenitis Suppurativa.")
    print("\nThe correct diagnosis is C. Hidradenitis Supportiva.")

# Execute the function to get the diagnosis
diagnose_skin_condition()