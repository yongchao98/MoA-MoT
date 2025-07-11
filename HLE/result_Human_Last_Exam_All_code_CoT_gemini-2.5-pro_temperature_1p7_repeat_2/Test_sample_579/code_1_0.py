def diagnose_skin_condition():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    # Step 1: Summarize the patient's key information.
    patient_profile = {
        "Locations": ["axillary folds", "inframammary folds", "inguinal regions"],
        "Lesions": ["large bullae", "erythematous plaques", "purulent nodules"],
        "Key Lesion": "purulent nodules",
        "Risk Factors": ["BMI of 39 (obesity)", "smoker"]
    }

    print("Starting diagnostic analysis based on the clinical vignette.")
    print("-" * 50)

    # Step 2: Analyze the location of the lesions.
    print("Analysis Step 1: Lesion Distribution")
    print(f"The lesions are found in the {', '.join(patient_profile['Locations'])}.")
    print("This distribution in skin folds (intertriginous areas) is characteristic of several conditions, including inverse psoriasis and Hidradenitis Suppurativa (HS).")
    print("-" * 50)

    # Step 3: Analyze the morphology of the lesions.
    print("Analysis Step 2: Lesion Morphology")
    print(f"The patient presents with a mix of lesions: {', '.join(patient_profile['Lesions'])}.")
    print(f"The most specific and diagnostically important lesion described is '{patient_profile['Key Lesion']}'.")
    print("Purulent (pus-filled) nodules are a hallmark of Hidradenitis Suppurativa.")
    print("While erythematous plaques can be seen in inverse psoriasis, and bullae can occur in various conditions, the purulent nodules are the strongest clue.")
    print("-" * 50)

    # Step 4: Analyze the patient's risk factors.
    print("Analysis Step 3: Patient Risk Factors")
    print(f"The patient has a {patient_profile['Risk Factors'][0]} and is a {patient_profile['Risk Factors'][1]}.")
    print("Both obesity and smoking are major, well-established risk factors for Hidradenitis Suppurativa.")
    print("-" * 50)

    # Step 5: Synthesize findings and conclude.
    print("Conclusion:")
    print("While elements of other conditions might be present (e.g., plaques resembling psoriasis), the combination of:")
    print("  1. Classic locations (axillae, inframammary, inguinal)")
    print("  2. The hallmark lesion (purulent nodules)")
    print("  3. Key patient risk factors (obesity, smoking)")
    print("makes Hidradenitis Suppurativa the most compelling diagnosis.")
    print("\nFinal Diagnosis Choice: C. Hidradenitis Supportiva")

# Run the diagnostic analysis
diagnose_skin_condition()