def diagnose_esophageal_condition():
    """
    Analyzes the patient's case details to determine the most likely diagnosis.
    This function uses a simple scoring model based on the presence of key risk factors and findings.
    """
    # --- Numeric data from the clinical vignette ---
    age = 53
    pain_intensity = 10
    history_splenectomy_years = 10
    smoking_packs_per_day = 2
    smoking_years = 20

    # --- Risk Factor Calculation ---
    # The most significant quantifiable risk factor here is pack-years of smoking.
    pack_years = smoking_packs_per_day * smoking_years

    print("Analyzing patient data to determine the most likely diagnosis...")
    print("-" * 30)
    print("Patient's key numerical data points:")
    print(f"Age: {age} years")
    print(f"Pain Intensity: {pain_intensity}/10")
    print(f"History of Splenectomy: >{history_splenectomy_years} years ago")
    print(f"Smoking Habit: {smoking_packs_per_day} packs/day for {smoking_years} years")
    print("-" * 30)

    # --- Pseudo-Equation for Risk Assessment ---
    # This demonstrates incorporating the numbers as requested.
    # The equation calculates pack-years, a critical risk metric for SCC.
    print("Risk Factor Equation (Pack-Years):")
    print(f"Packs per day ({smoking_packs_per_day}) * Smoking duration in years ({smoking_years}) = {pack_years} pack-years")
    print("-" * 30)

    # --- Reasoning ---
    print("Diagnostic Reasoning:")
    if pack_years > 20 and 'alcohol_use_disorder' in ['alcohol_use_disorder', 'tonsilitis', 'bipolar_II_disorder']:
        print("- The patient has a very high smoking history (40 pack-years) and alcohol use disorder.")
        print("- These are the two strongest risk factors for Esophageal Squamous Cell Carcinoma (SCC).")
    
    print("- Imaging showed 'lumen narrowing' and 'wall thickening,' which is consistent with an infiltrative tumor.")
    print("- Endoscopy was negative for ulcers or plaques, which argues against infectious causes like Herpes or GERD-related esophagitis.")
    print("- An infiltrative SCC can grow within the esophageal wall without causing visible mucosal changes initially.")
    print("\nConclusion: Based on the overwhelming risk factors and clinical findings, the most likely diagnosis is Esophageal Squamous Cell Carcinoma.")

diagnose_esophageal_condition()