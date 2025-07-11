def solve_diagnosis():
    """
    This function analyzes a clinical vignette to determine the most likely diagnosis
    by scoring each option based on key patient data.
    """

    # --- Patient's Key Clinical Features ---
    # We assign a value of 1 for presence of a feature, 0 for absence.
    # We use a higher value for major risk factors and key findings.
    risk_factors = {
        'heavy_smoking': 2,  # Major risk factor
        'alcohol_use_disorder': 2,  # Major risk factor
        'gerd_history': 0
    }

    findings = {
        'odynophagia': 1,  # Pain with swallowing
        'imaging_wall_thickening': 2,  # Key finding
        'endoscopy_normal_mucosa': 2   # Key finding (paradoxical)
    }

    # --- Initialize Scores for Each Diagnosis ---
    scores = {
        'A. Streptococcal esophagitis': 0,
        'B. Esophageal adenocarcinoma': 0,
        'C. Esophageal squamous cell carcinoma': 0,
        'D. GERD': 0,
        'E. Herpes esophagitis': 0,
    }

    print("Calculating diagnostic scores based on patient data...\n")

    # --- Score Calculation ---

    # A & E: Infectious Esophagitis (Strep & Herpes)
    # The normal endoscopy finding strongly argues against infectious causes.
    score_a = 0 - (findings['endoscopy_normal_mucosa'] * 3)
    scores['A. Streptococcal esophagitis'] = score_a
    print(f"Equation for Strep Esophagitis (A): 0 - (Normal Endoscopy * 3) = {score_a}")

    score_e = 0 - (findings['endoscopy_normal_mucosa'] * 3)
    scores['E. Herpes esophagitis'] = score_e
    print(f"Equation for Herpes Esophagitis (E): 0 - (Normal Endoscopy * 3) = {score_e}\n")


    # D: GERD
    # Normal endoscopy and wall thickening on imaging are not typical for GERD.
    score_d = 0 - findings['imaging_wall_thickening'] - findings['endoscopy_normal_mucosa']
    scores['D. GERD'] = score_d
    print(f"Equation for GERD (D): 0 - Wall Thickening ({findings['imaging_wall_thickening']}) - Normal Endoscopy ({findings['endoscopy_normal_mucosa']}) = {score_d}\n")


    # B: Esophageal Adenocarcinoma
    # Wrong risk factors (patient has smoking/alcohol, not GERD history).
    score_b = findings['imaging_wall_thickening'] - (risk_factors['heavy_smoking'] + risk_factors['alcohol_use_disorder'])
    scores['B. Esophageal adenocarcinoma'] = score_b
    print(f"Equation for Adenocarcinoma (B): Wall Thickening ({findings['imaging_wall_thickening']}) - Smoking ({risk_factors['heavy_smoking']}) - Alcohol ({risk_factors['alcohol_use_disorder']}) = {score_b}\n")


    # C: Esophageal Squamous Cell Carcinoma (SCC)
    # Strong positive points for risk factors (smoking, alcohol).
    # The combination of imaging (wall thickening) and a normal endoscopy (submucosal invasion) is a classic presentation for infiltrating SCC.
    score_c = (risk_factors['heavy_smoking'] * 2) + \
              (risk_factors['alcohol_use_disorder'] * 2) + \
              (findings['imaging_wall_thickening'] * findings['endoscopy_normal_mucosa'])
    scores['C. Esophageal squamous cell carcinoma'] = score_c
    print("Equation for Squamous Cell Carcinoma (C):")
    print(f"(Smoking * 2) + (Alcohol * 2) + (Wall Thickening * Normal Endoscopy) = "
          f"({risk_factors['heavy_smoking']} * 2) + ({risk_factors['alcohol_use_disorder']} * 2) + ({findings['imaging_wall_thickening']} * {findings['endoscopy_normal_mucosa']}) = {score_c}\n")

    # --- Determine the most likely diagnosis ---
    most_likely_diagnosis = max(scores, key=scores.get)
    print("--- Conclusion ---")
    print(f"The diagnosis with the highest score is: {most_likely_diagnosis}")


solve_diagnosis()
<<<C>>>