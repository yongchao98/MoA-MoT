def find_most_likely_diagnosis():
    """
    Analyzes clinical findings to determine the most likely diagnosis
    by scoring each option against the patient's symptoms.
    """
    # Patient's clinical findings from the vignette.
    # Spasticity is an upper motor neuron sign (stiffness), distinct from the
    # muscle weakness (flaccidity) typically seen in myositis.
    # Hypertrophic scarring can be a result of chronic skin inflammation or calcinosis.
    # Anti-Mi-2 is a specific antibody for dermatomyositis, but is often negative in the juvenile form.
    patient_findings = {
        "hypertrophic_scarring": True,
        "erythema": True,
        "spasticity": True,
        "negative_anti_mi2": True
    }

    # Knowledge base scoring symptoms against each possible diagnosis.
    # +1.0: Classic match, +0.5: Plausible association, 0.0: Unrelated, -1.0: Contradiction
    diagnosis_scores = {
        "A. Ectropion": {
            "name": "Ectropion",
            "hypertrophic_scarring": 0.0,
            "erythema": 0.0,
            "spasticity": 0.0,
            "negative_anti_mi2": 0.0
        },
        "B. McArdle disease": {
            "name": "McArdle disease",
            "hypertrophic_scarring": -1.0,
            "erythema": -1.0,
            "spasticity": -1.0,
            "negative_anti_mi2": 0.0
        },
        "C. Dermatomyositis": {
            "name": "Dermatomyositis (Juvenile form considered)",
            "hypertrophic_scarring": 0.5, # Plausible from calcinosis cutis, which can ulcerate and scar
            "erythema": 1.0,              # A classic feature (vasculitis-related rash)
            "spasticity": -1.0,           # A major contradiction (weakness is typical, but severe contractures can mimic spasticity)
            "negative_anti_mi2": 0.5      # Consistent with Juvenile Dermatomyositis
        },
        "D. McCune Albright syndrome": {
            "name": "McCune Albright syndrome",
            "hypertrophic_scarring": -1.0,
            "erythema": -1.0,
            "spasticity": 0.0,
            "negative_anti_mi2": 0.0
        },
        "E. Cataracts": {
            "name": "Cataracts",
            "hypertrophic_scarring": -1.0,
            "erythema": -1.0,
            "spasticity": -1.0,
            "negative_anti_mi2": 0.0
        }
    }

    print("Calculating scores for each possible diagnosis based on patient findings...\n")
    
    best_diagnosis = None
    max_score = -float('inf')

    for key, data in diagnosis_scores.items():
        score_h = data["hypertrophic_scarring"]
        score_e = data["erythema"]
        score_s = data["spasticity"]
        score_l = data["negative_anti_mi2"]
        
        total_score = score_h + score_e + score_s + score_l

        if total_score > max_score:
            max_score = total_score
            best_diagnosis = key

        # Output the calculation for each diagnosis as requested
        print(f"Diagnosis: {key}")
        print(f"Equation: {score_h} (scarring) + {score_e} (erythema) + {score_s} (spasticity) + {score_l} (lab) = {total_score:.1f}")
        print("-" * 30)

    print(f"\nThe diagnosis with the highest score is '{best_diagnosis}', making it the most likely choice.")

find_most_likely_diagnosis()