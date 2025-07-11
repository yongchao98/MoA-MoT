import collections

def diagnose_patient():
    """
    Analyzes patient symptoms to find the most likely diagnosis from a list of options.
    """
    # Step 1: Define the key patient findings from the clinical case.
    # The patient has a pelvic mass, hypertension, aniridia, and developmental delay (delayed speech).
    patient_findings = {"pelvic_mass", "elevated_blood_pressure", "aniridia", "developmental_delay"}

    # Step 2: Create a knowledge base of disease associations.
    disease_associations = {
        "A. Germ cell tumor": ["pelvic_mass"],
        "B. Astrocytoma": ["neurological_deficits"], # Primarily a brain tumor
        "C. Neuroblastoma": ["pelvic_mass", "elevated_blood_pressure"],
        "D. Nephroblastoma": ["pelvic_mass", "elevated_blood_pressure", "aniridia", "developmental_delay"],
        "E. Ewing sarcoma": ["pelvic_mass"] # Primarily a bone tumor
    }

    print("Analyzing patient findings against possible diagnoses...")
    print("--------------------------------------------------")
    print(f"Key patient findings: {list(patient_findings)}\n")

    scores = collections.defaultdict(int)

    # Step 3 & 4: Calculate a match score for each diagnosis.
    for diagnosis, associations in disease_associations.items():
        print(f"Evaluating: {diagnosis}")
        
        match_count = 0
        matching_symptoms = []

        # Count the number of patient findings explained by the disease
        for finding in patient_findings:
            if finding in associations:
                match_count += 1
                matching_symptoms.append(finding)
        
        base_score = match_count
        bonus_score = 0

        # Apply a special bonus for the WAGR syndrome association (Wilms Tumor/Nephroblastoma + Aniridia)
        if "Nephroblastoma" in diagnosis and "aniridia" in matching_symptoms:
            bonus_score = 5
            print("  - SPECIAL BONUS: Detected WAGR syndrome component (Nephroblastoma + aniridia). Adding +5.")

        final_score = base_score + bonus_score
        scores[diagnosis] = final_score

        # This section fulfills the requirement to "output each number in the final equation"
        print(f"  - Matching Symptoms: {len(matching_symptoms)} ({', '.join(matching_symptoms) if matching_symptoms else 'None'})")
        print(f"  - Score Calculation: {base_score} (matches) + {bonus_score} (bonus) = {final_score}")
        print("-" * 25)

    # Step 5: Determine the diagnosis with the highest score.
    best_diagnosis = max(scores, key=scores.get)
    
    print("\n---FINAL CONCLUSION---")
    print(f"The most likely diagnosis is '{best_diagnosis}' with a score of {scores[best_diagnosis]}.")
    print("\nREASONING:")
    print("Nephroblastoma (Wilms tumor) is the only diagnosis that accounts for the entire constellation of the patient's findings, including the pelvic mass, hypertension, and crucially, the aniridia and developmental delay which are classic components of WAGR syndrome (Wilms tumor, Aniridia, Genitourinary anomalies, and Retardation).")

diagnose_patient()
<<<D>>>