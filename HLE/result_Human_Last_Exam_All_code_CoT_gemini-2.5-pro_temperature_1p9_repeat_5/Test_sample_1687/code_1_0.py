import sys

def calculate_diagnosis_likelihood():
    """
    This script analyzes clinical findings from a patient case to determine the most likely diagnosis.
    It assigns points to each possible diagnosis based on the presence of key signs and symptoms.
    """
    # Key clinical findings extracted from the case vignette
    # A significant hemoglobin drop is defined here as > 5 g/dL (11.7 to 6.5)
    findings = {
        'procedure': 'difficult_colonoscopy',
        'no_polypectomy': True,
        'pain_location': 'LUQ',
        'referred_pain': 'left_shoulder',
        'hemodynamic_instability': True,
        'significant_hemoglobin_drop': True
    }

    # Define the possible diagnoses
    diagnoses = {
        "A. Colonic perforation": 0,
        "B. Lower GI bleeding": 0,
        "C. Splenic laceration": 0,
        "D. Postpolypectomy syndrome": 0,
    }

    print("Analyzing clinical data to determine the most likely diagnosis...")
    print("-" * 30)

    # --- Scoring Logic ---

    # Score for Splenic Laceration
    score_c = 0
    # Splenic injury is strongly associated with LUQ pain
    if findings['pain_location'] == 'LUQ':
        score_c += 3
    # Left shoulder pain (Kehr's sign) is highly specific for diaphragmatic irritation from splenic bleed
    if findings['referred_pain'] == 'left_shoulder':
        score_c += 3
    # Splenic artery is high flow, laceration causes massive bleeding
    if findings['significant_hemoglobin_drop']:
        score_c += 2
    # Massive blood loss leads to shock
    if findings['hemodynamic_instability']:
        score_c += 2
    diagnoses["C. Splenic laceration"] = score_c

    # Score for Colonic Perforation
    score_a = 0
    # Perforation can cause pain, but not specifically LUQ/shoulder combo
    if findings['pain_location'] == 'LUQ':
        score_a += 1
    # Perforation can cause some bleeding, but rarely this massive
    if findings['significant_hemoglobin_drop']:
        score_a += 1
    # Can lead to shock (septic or hypovolemic)
    if findings['hemodynamic_instability']:
        score_a += 1
    diagnoses["A. Colonic perforation"] = score_a

    # Score for Lower GI bleeding
    score_b = 0
    # This is a symptom, so it matches the bleeding finding, but it's not a specific diagnosis.
    # The pain location is also atypical for a "lower" source.
    if findings['significant_hemoglobin_drop']:
        score_b += 2
    if findings['hemodynamic_instability']:
        score_b += 1
    diagnoses["B. Lower GI bleeding"] = score_b

    # Score for Postpolypectomy syndrome
    score_d = 0
    # This diagnosis is ruled out by the history
    if findings['no_polypectomy']:
        score_d = -5 # Assign a negative score as it's explicitly contradicted
    diagnoses["D. Postpolypectomy syndrome"] = score_d

    # --- Print Results ---
    print("Likelihood scores based on findings:")
    for diagnosis, score in diagnoses.items():
        print(f"{diagnosis}: {score} points")

    most_likely_diagnosis = max(diagnoses, key=diagnoses.get)
    
    print("-" * 30)
    print(f"Conclusion: The most likely diagnosis is '{most_likely_diagnosis}'.")
    print("This diagnosis perfectly aligns the LUQ pain, referred left shoulder pain (Kehr's sign), and evidence of massive intra-abdominal hemorrhage following a difficult colonoscopy.")
    
    # Fulfilling the request to show the final equation for the highest-scoring diagnosis
    print("\nFinal Score Calculation for the Winning Diagnosis:")
    c_luq = 3 if findings['pain_location'] == 'LUQ' else 0
    c_shoulder = 3 if findings['referred_pain'] == 'left_shoulder' else 0
    c_hgb = 2 if findings['significant_hemoglobin_drop'] else 0
    c_instability = 2 if findings['hemodynamic_instability'] else 0
    
    # We output each number in the final equation.
    print(f"'{most_likely_diagnosis}' Score = {c_luq} (LUQ Pain) + {c_shoulder} (Left Shoulder Pain) + {c_hgb} (Hgb Drop) + {c_instability} (Instability)")
    print(f"Total Score = {diagnoses[most_likely_diagnosis]}")


if __name__ == "__main__":
    calculate_diagnosis_likelihood()
    sys.stdout.flush() # Ensure all print statements are shown before the final answer tag.
    print("<<<C>>>")