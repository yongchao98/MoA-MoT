def diagnose_patient():
    """
    Analyzes patient features against potential diagnoses to find the best fit.
    This is a simplified model for demonstrating clinical reasoning.
    """
    # Patient's key clinical features
    patient_features = [
        "Short Stature",
        "Ovarian Dysgenesis",
        "Exercise Intolerance",
        "Hypertension",
        "Normal Karyotype"
    ]

    # Database of diagnoses and their associated features
    # Scoring: 1 = strongly explains, 0.5 = sometimes explains, 0 = does not explain, -100 = contradicts
    diagnoses_features = {
        "Turner Syndrome (Classic 45,X)": {
            "Short Stature": 1,
            "Ovarian Dysgenesis": 1,
            "Exercise Intolerance": 1,
            "Hypertension": 1,
            "Normal Karyotype": -100  # Ruled out by this key finding
        },
        "FMR1 Premutation (FXPOI)": {
            "Short Stature": 0.5,
            "Ovarian Dysgenesis": 1,
            "Exercise Intolerance": 0,
            "Hypertension": 0,
            "Normal Karyotype": 1
        },
        "Mitochondrial DNA Mutation": {
            "Short Stature": 1,
            "Ovarian Dysgenesis": 1,
            "Exercise Intolerance": 1,
            "Hypertension": 1,
            "Normal Karyotype": 1
        },
        "Isolated Gonadal Dysgenesis (e.g., NR5A1 mutation)": {
            "Short Stature": 0.5,
            "Ovarian Dysgenesis": 1,
            "Exercise Intolerance": 0,
            "Hypertension": 0,
            "Normal Karyotype": 1
        }
    }

    scores = {}
    print("Evaluating potential diagnoses based on patient features...\n")

    # Calculate score for each diagnosis
    for diagnosis, features in diagnoses_features.items():
        total_score = 0
        equation_parts = []
        for feature in patient_features:
            score = features.get(feature, 0)
            total_score += score
            equation_parts.append(str(score))
        
        scores[diagnosis] = total_score
        equation_str = " + ".join(equation_parts)
        print(f"Diagnosis: {diagnosis}")
        print(f"Scoring based on features ({', '.join(patient_features)}):")
        print(f"Equation: {equation_str} = {total_score}\n")

    # Find the best diagnosis
    # We filter out any diagnosis with a score below zero, as it contradicts a key finding.
    valid_scores = {dx: sc for dx, sc in scores.items() if sc > 0}
    if not valid_scores:
        best_diagnosis = "No suitable diagnosis found."
    else:
        best_diagnosis = max(valid_scores, key=valid_scores.get)

    print("---------------------------------------------------------")
    print(f"Conclusion: The most likely diagnosis is '{best_diagnosis}'.")
    print("This is because it provides the most comprehensive explanation for the patient's multi-systemic symptoms (involving growth, ovarian function, and muscle function) in the context of a normal karyotype.")

if __name__ == "__main__":
    diagnose_patient()