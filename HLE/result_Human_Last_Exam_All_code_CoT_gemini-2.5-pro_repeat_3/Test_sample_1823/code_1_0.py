def diagnose():
    """
    This function analyzes clinical findings to determine the most likely diagnosis.
    It uses a scoring system to illustrate the reasoning process, where a positive
    score indicates a feature is consistent with a diagnosis, and a negative score
    indicates a feature is inconsistent.
    """
    # Clinical Findings from the Vignette
    # 1. Hypertrophic scarring
    # 2. Erythema
    # 3. Spasticity
    # 4. Negative anti-Mi-2

    print("Evaluating potential diagnoses based on the patient's symptoms...\n")

    # --- Diagnosis C: Dermatomyositis ---
    print("Analysis for C. Dermatomyositis:")
    score_c = 0
    
    # Score for 'Erythema' (skin inflammation is a key feature)
    score_c_erythema = 1
    score_c += score_c_erythema
    print(f"Finding 'Erythema' is consistent: +{score_c_erythema} point")

    # Score for 'Hypertrophic scarring' (can result from calcinosis cutis, a known complication)
    score_c_scarring = 1
    score_c += score_c_scarring
    print(f"Finding 'Hypertrophic scarring' is consistent: +{score_c_scarring} point")

    # Score for 'Spasticity' (interpreted as severe stiffness/contractures from myositis)
    score_c_spasticity = 1
    score_c += score_c_spasticity
    print(f"Finding 'Spasticity' is a plausible fit: +{score_c_spasticity} point")
    
    # Score for 'Negative anti-Mi-2' (common in Juvenile Dermatomyositis, so it doesn't rule it out)
    score_c_labs = 0
    score_c += score_c_labs
    print(f"Finding 'Negative anti-Mi-2' is neutral: +{score_c_labs} points")

    print(f"Final Equation and Score for Dermatomyositis: {score_c_erythema} + {score_c_scarring} + {score_c_spasticity} + {score_c_labs} = {score_c}\n")

    # --- Diagnosis B: McArdle disease ---
    print("Analysis for B. McArdle disease:")
    score_b = 0

    # Score for 'Erythema' (not a feature)
    score_b_erythema = -1
    score_b += score_b_erythema
    print(f"Finding 'Erythema' is inconsistent: {score_b_erythema} point")

    # Score for 'Hypertrophic scarring' (not a feature)
    score_b_scarring = -1
    score_b += score_b_scarring
    print(f"Finding 'Hypertrophic scarring' is inconsistent: {score_b_scarring} point")

    # Score for 'Spasticity' (inconsistent; weakness/fatigue is typical)
    score_b_spasticity = -1
    score_b += score_b_spasticity
    print(f"Finding 'Spasticity' is inconsistent: {score_b_spasticity} point")

    print(f"Final Equation and Score for McArdle disease: {score_b_erythema} + {score_b_scarring} + {score_b_spasticity} = {score_b}\n")
    
    print("Conclusion: With the highest score, Dermatomyositis is the most likely diagnosis as it accounts for the full constellation of symptoms.")

diagnose()