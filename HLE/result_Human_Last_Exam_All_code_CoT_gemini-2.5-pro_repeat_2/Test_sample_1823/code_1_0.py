def find_diagnosis():
    """
    This function simulates a diagnostic scoring process to determine the most likely diagnosis
    based on the patient's clinical features.
    """
    patient_features = {
        "age_infant": 1,
        "hypertrophic_scarring": 1,
        "erythema": 1,
        "spasticity": 1,
        "negative_anti_mi2": 1
    }

    # Scoring matrix: How well each feature fits a diagnosis.
    # Higher score means better fit.
    diagnoses = {
        "A. Ectropion": {
            "age_infant": 0, "hypertrophic_scarring": 0, "erythema": 0, "spasticity": -1, "negative_anti_mi2": 0
        },
        "B. McArdle disease": {
            "age_infant": -1, "hypertrophic_scarring": -1, "erythema": 0, "spasticity": -2, "negative_anti_mi2": 0
        },
        "C. Dermatomyositis": {
            "age_infant": 1, "hypertrophic_scarring": 0.5, "erythema": 2, "spasticity": 0.5, "negative_anti_mi2": 1
        },
        "D. McCune Albright syndrome": {
            "age_infant": 1, "hypertrophic_scarring": -1, "erythema": -1, "spasticity": -1, "negative_anti_mi2": 0
        },
        "E. Cataracts": {
            "age_infant": 0, "hypertrophic_scarring": -1, "erythema": -1, "spasticity": -1, "negative_anti_mi2": 0
        }
    }

    best_diagnosis = None
    max_score = -float('inf')
    final_equation_components = {}

    print("Calculating diagnostic scores...\n")

    for diagnosis, scores in diagnoses.items():
        total_score = sum(scores.values())
        print(f"Diagnosis: {diagnosis}, Score: {total_score}")
        if total_score > max_score:
            max_score = total_score
            best_diagnosis = diagnosis
            final_equation_components = scores

    print("\n--- Conclusion ---")
    print(f"The most likely diagnosis is: {best_diagnosis}")
    print("This is based on the following scoring calculation:")

    # Retrieve the scores for the best diagnosis
    s = final_equation_components
    equation_str = (
        f"{s['age_infant']} (for infant age) + "
        f"{s['hypertrophic_scarring']} (for hypertrophic scarring) + "
        f"{s['erythema']} (for erythema) + "
        f"{s['spasticity']} (for spasticity) + "
        f"{s['negative_anti_mi2']} (for negative anti-Mi-2) = {max_score}"
    )

    print(f"\nFinal Equation for {best_diagnosis.split('. ')[1]}:\n{equation_str}")

find_diagnosis()
<<<C>>>