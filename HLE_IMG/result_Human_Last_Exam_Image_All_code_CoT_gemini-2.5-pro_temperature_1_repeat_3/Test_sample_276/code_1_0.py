def calculate_diagnostic_score():
    """
    This function calculates a simplified score for potential diagnoses
    based on the patient's clinical features.
    Features are weighted based on their diagnostic significance.
    """

    # Clinical features from the case
    features = {
        "Ileocecal_Involvement": 1,
        "Chronic_Course (>4wks)": 1,
        "Constitutional_Symptoms (Fever/Wt_Loss)": 1,
        "Extra-intestinal_Manifestations (Uveitis/Arthritis)": 2, # Highly specific for IBD
        "Acute_Pain/Inflammation": 1,
        "Advanced_Age": 1
    }

    # Scoring for top differential diagnoses
    crohns_score = (features["Ileocecal_Involvement"] +
                   features["Chronic_Course (>4wks)"] +
                   features["Constitutional_Symptoms (Fever/Wt_Loss)"] +
                   features["Extra-intestinal_Manifestations (Uveitis/Arthritis)"] +
                   features["Acute_Pain/Inflammation"] +
                   features["Advanced_Age"])

    tb_score = (features["Ileocecal_Involvement"] +
                features["Chronic_Course (>4wks)"] +
                features["Constitutional_Symptoms (Fever/Wt_Loss)"] +
                # TB doesn't typically cause this combination of EIMs
                0 +
                features["Acute_Pain/Inflammation"] +
                features["Advanced_Age"])

    malignancy_score = (features["Ileocecal_Involvement"] +
                        features["Chronic_Course (>4wks)"] +
                        features["Constitutional_Symptoms (Fever/Wt_Loss)"] +
                        # Malignancy does not cause these EIMs
                        0 +
                        features["Acute_Pain/Inflammation"] +
                        features["Advanced_Age"])

    print("Diagnostic Reasoning Score Calculation:")
    print("-" * 40)
    print(f"A. Crohn's Disease Score: 1 + 1 + 1 + 2 + 1 + 1 = {crohns_score}")
    print(f"C. Ileocecal TB Score: 1 + 1 + 1 + 0 + 1 + 1 = {tb_score}")
    print(f"K. GI Malignancy Score: 1 + 1 + 1 + 0 + 1 + 1 = {malignancy_score}")
    print("-" * 40)
    print("Conclusion: Crohn's Disease has the highest score as it best explains the complete clinical picture, including the crucial history of uveitis and arthritis.")

calculate_diagnostic_score()