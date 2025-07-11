def diagnose_case():
    """
    Analyzes clinical findings from a patient case to determine the most likely diagnosis.
    The function uses a simple scoring system to weigh the evidence for each potential diagnosis.
    """
    # Key clinical findings extracted from the case description
    findings = {
        'procedure': 'Difficult colonoscopy',
        'polypectomy_performed': False,
        'pain_location': 'Left Upper Quadrant and epigastrium',
        'referred_pain': 'Left-shoulder discomfort (Kehr\'s sign)',
        'initial_hemoglobin': 11.7,
        'later_hemoglobin': 6.5,
        'hemodynamic_status': 'Tachycardia and hypotension (instability)'
    }

    # Initialize scores for each diagnosis
    scores = {
        'A. Colonic perforation': 0,
        'B. Lower GI bleeding': 0,
        'C. Splenic laceration': 0,
        'D. Postpolypectomy syndrome': 0,
    }

    print("Analyzing clinical findings to determine the most likely diagnosis...")
    print("-" * 30)

    # --- Scoring Logic ---

    # Score for A. Colonic perforation
    # A perforation can cause instability, and a high tear could cause LUQ pain.
    scores['A. Colonic perforation'] += 2  # Consistent with difficult colonoscopy
    scores['A. Colonic perforation'] += 1  # Pain and instability are possible, but not the most classic signs

    # Score for B. Lower GI bleeding
    # Inconsistent pain location and lack of a polypectomy make this less likely.
    scores['B. Lower GI bleeding'] -= 2 # Pain location is atypical (upper, not lower abdomen)
    if not findings['polypectomy_performed']:
        scores['B. Lower GI bleeding'] -= 3 # Major cause (polypectomy) is absent

    # Score for C. Splenic laceration
    # This diagnosis strongly aligns with the specific combination of findings.
    hemoglobin_drop = findings['initial_hemoglobin'] - findings['later_hemoglobin']
    
    scores['C. Splenic laceration'] += 2 # Consistent with difficult colonoscopy
    scores['C. Splenic laceration'] += 3 # Classic pain location (LUQ)
    scores['C. Splenic laceration'] += 5 # Highly specific referred pain (Kehr's sign to left shoulder)
    scores['C. Splenic laceration'] += 5 # Massive hemorrhage (Hemoglobin drop of 5.2) and instability

    # Score for D. Postpolypectomy syndrome
    # This is ruled out by the patient's history.
    if not findings['polypectomy_performed']:
        scores['D. Postpolypectomy syndrome'] -= 10 # Ruled out as no polypectomy was done

    # --- Print Results ---
    print("Scoring each diagnosis based on the evidence:")
    for diagnosis, score in scores.items():
        # This loop demonstrates printing numbers from our 'equation' (scoring system)
        print(f"Final Score for {diagnosis}: {score}")

    print("-" * 30)

    # Determine the most likely diagnosis
    most_likely_diagnosis = max(scores, key=scores.get)
    
    # Extract just the letter for the final answer format
    final_answer_letter = most_likely_diagnosis.split('.')[0]

    print(f"The diagnosis with the highest score is: {most_likely_diagnosis}")
    print(f"\nThe patient's presentation of left upper quadrant pain, referred left shoulder pain (Kehr's sign), and rapid hemorrhagic shock after a difficult colonoscopy is classic for a splenic injury.")
    print(f"\n<<<C>>>")

if __name__ == "__main__":
    diagnose_case()