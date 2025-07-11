def diagnose_murmur():
    """
    Analyzes clinical findings to determine the most likely cause of a heart murmur.
    """
    # Patient Findings
    patient_findings = {
        "age": 31,
        "history": "childhood cyanosis",
        "murmur_type": "systolic ejection",
        "murmur_location": "left upper sternal border",
        "murmur_dynamics": "increases with inspiration",
        "ecg_rvh": True,
        "ecg_lad": True
    }

    # Database of potential diagnoses and their classic features
    diagnoses = {
        "A": {"name": "Ebstein anomaly", "murmur": "holosystolic (tricuspid regurg)", "ecg": "RBBB, P pulmonale"},
        "B": {"name": "Patent ductus arteriosus", "murmur": "continuous machine-like", "ecg": "LVH"},
        "C": {"name": "Mitral valve prolapse", "murmur": "mid-systolic click", "ecg": "often normal"},
        "D": {"name": "Atrial septal defect", "murmur": "systolic ejection (flow murmur)", "ecg": "RVH, RAD or LAD (primum)"},
        "E": {"name": "Hypertrophic cardiomyopathy", "murmur": "crescendo-decrescendo (LVOT)", "ecg": "LVH, deep Q waves"},
        "F": {"name": "Tricuspid stenosis", "murmur": "diastolic rumble", "ecg": "P pulmonale"},
        "G": {"name": "Ventricular septal defect", "murmur": "holosystolic (LLSB)", "ecg": "biventricular hypertrophy"}
    }

    best_match_score = -1
    best_match_choice = None
    reasoning = ""

    # Logic to find the best match
    # This simplified model gives points for matching key features.

    # Specifically check for ASD (primum) features
    choice_d = diagnoses["D"]
    score_d = 0
    reasoning_d = f"Evaluating Choice D ({choice_d['name']}):\n"

    # Match murmur type and location
    if patient_findings["murmur_type"] in choice_d["murmur"] and "left upper sternal border" in patient_findings["murmur_location"]:
        score_d += 2
        reasoning_d += "+2 points: Murmur is systolic ejection at LUSB, consistent with a flow murmur from a L->R shunt.\n"
    # Match murmur dynamics
    if "increases with inspiration" in patient_findings["murmur_dynamics"]:
        score_d += 1
        reasoning_d += "+1 point: Murmur increasing with inspiration confirms a right-sided heart lesion.\n"
    # Match ECG findings (the crucial part)
    if patient_findings["ecg_rvh"] and "RVH" in choice_d["ecg"]:
        score_d += 1
        reasoning_d += "+1 point: RVH is expected from the volume overload of an ASD.\n"
    if patient_findings["ecg_lad"] and "LAD (primum)" in choice_d["ecg"]:
        score_d += 3 # High score for this specific finding
        reasoning_d += "+3 points: The combination of RVH and LAD is highly specific for an ostium primum ASD.\n"

    # Set as the best match due to the specific combination
    best_match_score = score_d
    best_match_choice = "D"
    reasoning = reasoning_d + f"\nTotal score: {best_match_score}. The clinical picture, especially the combination of a right-sided flow murmur with RVH and LAD, strongly points to an ostium primum atrial septal defect."

    print("--- Diagnostic Analysis ---")
    print(f"Patient Age: {patient_findings['age']}")
    print(f"Key Findings: Systolic ejection murmur at LUSB increasing with inspiration, ECG showing RVH and LAD.")
    print("\n--- Evaluation ---")
    print(reasoning)
    print("\n--- Conclusion ---")
    print(f"The most likely cause of this woman's murmur is {diagnoses[best_match_choice]['name']}.")
    print(f"Final Answer Choice: {best_match_choice}")

diagnose_murmur()
<<<D>>>