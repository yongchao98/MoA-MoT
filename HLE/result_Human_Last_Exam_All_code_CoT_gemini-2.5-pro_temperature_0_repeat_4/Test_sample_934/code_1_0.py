def diagnose_esophageal_condition():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    This function simulates a diagnostic reasoning process by scoring potential diagnoses
    based on the patient's signs, symptoms, and risk factors.
    """

    # Patient Data from the vignette
    patient_profile = {
        "age": 53,
        "symptoms": ["10/10 substernal chest pain", "pain with swallowing (odynophagia)"],
        "risk_factors": {
            "alcohol_use_disorder": True,
            "heavy_smoker": True, # 2 packs/day for 20 years
        },
        "labs": {
            "c_reactive_protein": "elevated",
            "leukocytosis": True,
        },
        "imaging": {
            "description": "esophageal lumen narrowing, wall thickening",
        },
        "endoscopy": {
            "description": "no signs of erythema, ulcers, plaques, or strictures",
        }
    }

    # Initialize scores for each diagnosis
    scores = {
        "A. Streptococcal esophagitis": 0,
        "B. Esophageal adenocarcinoma": 0,
        "C. Esophageal squamous cell carcinoma": 0,
        "D. GERD": 0,
        "E. Herpes esophagitis": 0,
    }

    print("Analyzing clinical findings...\n")

    # --- Scoring Logic ---

    # 1. Risk Factors: Heavy smoking and alcohol are major risks for SCC.
    if patient_profile["risk_factors"]["alcohol_use_disorder"]:
        scores["C. Esophageal squamous cell carcinoma"] += 3
        scores["B. Esophageal adenocarcinoma"] += 1
    if patient_profile["risk_factors"]["heavy_smoker"]:
        scores["C. Esophageal squamous cell carcinoma"] += 3
        scores["B. Esophageal adenocarcinoma"] += 1
    print("Patient has major risk factors for Squamous Cell Carcinoma (heavy smoking and alcohol use).")
    print(f"Current Scores -> SCC: {scores['C. Esophageal squamous cell carcinoma']}, Adenocarcinoma: {scores['B. Esophageal adenocarcinoma']}\n")


    # 2. Imaging Findings: Wall thickening/narrowing suggests an infiltrative process.
    scores["C. Esophageal squamous cell carcinoma"] += 3
    scores["B. Esophageal adenocarcinoma"] += 2
    print("Imaging shows wall thickening and narrowing, consistent with an infiltrative cancer.")
    print(f"Current Scores -> SCC: {scores['C. Esophageal squamous cell carcinoma']}, Adenocarcinoma: {scores['B. Esophageal adenocarcinoma']}\n")

    # 3. Endoscopy Findings: This is the key differentiator.
    # A normal mucosa rules out diagnoses that MUST have visible lesions.
    # It strongly suggests a submucosal/infiltrative process, classic for some SCCs.
    scores["A. Streptococcal esophagitis"] -= 5 # Would have exudates/ulcers
    scores["E. Herpes esophagitis"] -= 5 # Would have classic "volcano" ulcers
    scores["D. GERD"] -= 3 # Severe GERD causing these symptoms/imaging would show esophagitis
    scores["B. Esophageal adenocarcinoma"] -= 4 # Almost always presents as a visible mass/ulcer
    scores["C. Esophageal squamous cell carcinoma"] += 3 # This diagnosis uniquely explains normal mucosa with abnormal imaging (infiltrative growth)
    print("Crucial Finding: Endoscopy is normal. This rules out conditions with typical mucosal lesions (like Herpes or severe GERD).")
    print("However, an infiltrative squamous cell carcinoma can grow within the esophageal wall, causing thickening seen on imaging while the inner lining (mucosa) appears normal on endoscopy.")
    print("This finding strongly points towards an infiltrative cancer.\n")

    # --- Final Tally ---
    print("--- Final Score Calculation ---")
    for diagnosis, score in scores.items():
        print(f"Diagnosis: {diagnosis}, Final Score: {score}")

    # Determine the most likely diagnosis
    most_likely_diagnosis = max(scores, key=scores.get)

    print("\n--- Conclusion ---")
    print(f"The diagnosis that best fits all the clinical evidence, especially the combination of major risk factors and the discrepancy between imaging and endoscopy, is: {most_likely_diagnosis}")

    # Extract the letter for the final answer format
    final_answer_letter = most_likely_diagnosis.split('.')[0]
    return final_answer_letter

# Run the analysis and print the final answer
final_answer = diagnose_esophageal_condition()
print(f"\n<<<C>>>")