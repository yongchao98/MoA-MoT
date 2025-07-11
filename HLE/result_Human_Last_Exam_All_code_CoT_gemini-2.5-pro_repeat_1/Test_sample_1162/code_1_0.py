def diagnose_patient():
    """
    Analyzes patient symptoms to determine the most likely diagnosis from a list of options.
    """
    patient_findings = {
        "age": 2,
        "presentation": "lethargy",
        "vitals": "hypertension",
        "growth": "10th percentile (delayed)",
        "exam": ["conjunctival pallor (anemia)", "aniridia", "delayed speech"],
        "imaging": "pelvic mass"
    }

    # Simplified knowledge base of associated findings for each diagnosis
    # Note: Aniridia is a very specific finding.
    disease_profiles = {
        "A. Germ cell tumor": ["pelvic mass"],
        "B. Astrocytoma": ["brain tumor", "neurological deficits"], # Mismatch: Location
        "C. Neuroblastoma": ["abdominal/pelvic mass", "hypertension", "common in young children"],
        "D. Nephroblastoma (Wilms Tumor)": ["abdominal/pelvic mass", "hypertension", "common in young children", "anemia", "associated with WAGR syndrome (Wilms, Aniridia, Genitourinary, Retardation)"],
        "E. Ewing sarcoma": ["bone or soft tissue mass", "pelvic mass"]
    }

    print("Analyzing Patient Case:")
    print(f"- Age: {patient_findings['age']} years old")
    print(f"- Key Findings: {', '.join(patient_findings['exam'])}")
    print(f"- Vitals: {patient_findings['vitals']}")
    print(f"- Imaging: {patient_findings['imaging']}\n")

    best_match = None
    highest_score = -1
    match_reasoning = ""

    # Scoring the match
    for diagnosis, features in disease_profiles.items():
        score = 0
        current_reasoning = []
        
        # Check for matching key features
        if "pelvic mass" in patient_findings["imaging"] and ("pelvic mass" in features or "abdominal/pelvic mass" in features):
            score += 1
            current_reasoning.append("explains the pelvic mass")
        
        if patient_findings["vitals"] == "hypertension" and "hypertension" in features:
            score += 1
            current_reasoning.append("explains hypertension")

        if "conjunctival pallor (anemia)" in patient_findings["exam"] and "anemia" in features:
            score += 1
            current_reasoning.append("explains anemia (pallor)")

        # The key differentiating feature: Aniridia
        if "aniridia" in patient_findings["exam"] and "associated with WAGR syndrome (Wilms, Aniridia, Genitourinary, Retardation)" in features:
            score += 5 # Heavily weight this pathognomonic finding
            current_reasoning.append("EXPLAINS THE HIGHLY SPECIFIC FINDING OF ANIRIDIA")
        
        if "delayed speech" in patient_findings["exam"] and "associated with WAGR syndrome (Wilms, Aniridia, Genitourinary, Retardation)" in features:
            score += 2 # WAGR includes developmental delay (Retardation)
            current_reasoning.append("explains developmental delay (delayed speech)")


        if score > highest_score:
            highest_score = score
            best_match = diagnosis
            match_reasoning = f"'{best_match}' is the most likely diagnosis because it {', and '.join(current_reasoning)}."

    print("Diagnostic Reasoning:")
    print(match_reasoning)
    print("\nThe constellation of aniridia, a pelvic/abdominal mass (Wilms tumor), hypertension, and developmental delay points strongly to WAGR syndrome.")
    
    final_answer_choice = best_match.split('.')[0]
    print(f"\nFinal Answer Choice: {final_answer_choice}")

diagnose_patient()
<<<D>>>