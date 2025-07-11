def diagnose_skin_condition():
    """
    Analyzes clinical findings using a weighted scoring model to determine
    the most likely dermatological diagnosis.
    """
    print("Starting diagnostic analysis based on clinical findings...")
    print("A scoring model will be used to quantify the likelihood of Hidradenitis Suppurativa (HS).\n")

    # --- Patient Data & Risk Factors ---
    # BMI of 39 is class II obesity, a major risk factor for HS.
    bmi = 39
    # Patient is a smoker, another major risk factor for HS.
    is_smoker = True
    # Multiple classic locations for HS are involved.
    involved_locations = ["axillary folds", "inframammary folds", "inguinal regions"]
    # The presence of purulent nodules is a hallmark of HS.
    has_purulent_nodules = True

    # --- Scoring System ---
    hs_score = 0
    equation_components = []

    print("Constructing the diagnostic score equation:")

    # Score for Obesity (BMI > 30)
    if bmi > 30:
        score = 3
        hs_score += score
        equation_components.append(str(score))
        print(f"Finding: Obesity (BMI={bmi}). This is a strong risk factor. Score contribution: {score}")

    # Score for Smoking
    if is_smoker:
        score = 3
        hs_score += score
        equation_components.append(str(score))
        print(f"Finding: Smoking History. This is a strong risk factor. Score contribution: {score}")

    # Score for Location Involvement
    num_locations = len(involved_locations)
    if num_locations >= 2:
        score = 4
        hs_score += score
        equation_components.append(str(score))
        print(f"Finding: Involvement of {num_locations} characteristic locations. Score contribution: {score}")

    # Score for Hallmark Lesions (Purulent Nodules)
    if has_purulent_nodules:
        score = 5
        hs_score += score
        equation_components.append(str(score))
        print(f"Finding: Presence of purulent nodules. This is a hallmark sign. Score contribution: {score}")

    # --- Final Calculation and Conclusion ---
    print("\n----------------------------------------")
    final_equation_str = " + ".join(equation_components)
    print(f"The final diagnostic equation based on the findings is:")
    print(f"{final_equation_str} = {hs_score}")
    print("----------------------------------------\n")

    print("Conclusion:")
    print("The total score is significantly high, strongly indicating a diagnosis of Hidradenitis Suppurativa (HS).")
    print("This condition is characterized by chronic inflammation in apocrine gland-rich areas (axillae, inguinal, inframammary), and is strongly associated with obesity and smoking. The presence of purulent nodules is a key diagnostic feature.")
    print("\nCompared to other options:")
    print("- Psoriasis & Atopic Dermatitis typically do not present with purulent nodules.")
    print("- Allergic Contact Dermatitis is unlikely to appear in three separate sites with such varied lesions.")
    print("- Malignant Intertrigo is much rarer and less likely to present with this specific combination of findings.")

# Execute the diagnostic analysis
diagnose_skin_condition()