def analyze_clinical_case():
    """
    This script analyzes the clinical findings from the case to determine the most likely diagnosis.
    It uses a simple scoring system to weigh the evidence.
    """
    # Define diagnoses and initialize their scores to 0
    diagnoses = {
        "A. Colonic perforation": 0,
        "B. Lower GI bleeding": 0,
        "C. Splenic laceration": 0,
        "D. Postpolypectomy syndrome": 0,
    }

    # --- Evidence from the Case ---

    # 1. Procedure type: No polypectomy was performed.
    # This virtually rules out Postpolypectomy syndrome.
    diagnoses["D. Postpolypectomy syndrome"] -= 10

    # 2. Pain Location: Left upper quadrant (LUQ) and left-shoulder discomfort.
    # This is a classic presentation for splenic injury (Kehr's sign).
    diagnoses["C. Splenic laceration"] += 5
    diagnoses["A. Colonic perforation"] += 1  # Less likely location
    diagnoses["B. Lower GI bleeding"] += 0 # Atypical presentation

    # 3. Hemoglobin Drop: Decreased from 11.7 g/dL to 6.5 g/dL.
    # This indicates massive hemorrhage, a hallmark of splenic rupture.
    print("Initial Hemoglobin: 11.7")
    print("Post-bleed Hemoglobin: 6.5")
    print("Change in Hemoglobin:", round(11.7 - 6.5, 1))
    diagnoses["C. Splenic laceration"] += 5 # Highly indicative
    diagnoses["A. Colonic perforation"] += 2 # Possible, but less common to be this severe
    diagnoses["B. Lower GI bleeding"] += 2 # Possible, but presentation is wrong

    # 4. Hemodynamic status: Tachycardia (Heart Rate 105) and hypotension.
    # These are signs of hemorrhagic shock.
    print("Heart Rate:", 105)
    print("Blood Pressure: 100/65")
    diagnoses["C. Splenic laceration"] += 3
    diagnoses["A. Colonic perforation"] += 1
    diagnoses["B. Lower GI bleeding"] += 1

    # --- Conclusion ---
    print("\n--- Scoring Results ---")
    for diagnosis, score in diagnoses.items():
        print(f"'{diagnosis}': Score = {score}")

    most_likely_diagnosis = max(diagnoses, key=diagnoses.get)
    print("\nBased on the analysis, the most likely diagnosis is:")
    print(most_likely_diagnosis)

analyze_clinical_case()
<<<C>>>