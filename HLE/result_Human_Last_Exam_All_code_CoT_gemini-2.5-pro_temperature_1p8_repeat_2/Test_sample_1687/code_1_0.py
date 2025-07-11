import pandas as pd

def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring each option against the key patient findings.
    """

    # Step 1: Define the key clinical findings from the case.
    patient_findings = {
        'Procedure': 'Difficult colonoscopy',
        'Polypectomy Performed': False,
        'Symptom Onset': '12 hours post-procedure',
        'Pain Location': 'Left Upper Quadrant (LUQ) / Epigastrium',
        'Referred Pain': 'Left-shoulder discomfort (Kehr\'s sign)',
        'Hemodynamics': 'Progressed to tachycardia and hypotension (shock)',
        'Hemoglobin Drop': 'Significant and acute (11.7 -> 6.5 g/dL)',
        'Physical Exam': 'Worsening LUQ tenderness, guarding, peritoneal signs'
    }

    # Step 2: Define the characteristics of each potential diagnosis.
    diagnoses_features = {
        'A. Colonic perforation': {
            'Consistent': ['Difficult colonoscopy', 'Peritoneal signs'],
            'Inconsistent': ['LUQ/epigastric + left shoulder pain (less typical)', 'Massive bleed from perforation is less common than from solid organ injury'],
            'Score': 0
        },
        'B. Lower GI bleeding': {
            'Consistent': ['Hemoglobin Drop', 'Can be a complication of colonoscopy'],
            'Inconsistent': ['Significant abdominal pain (many LGI bleeds are painless)', 'Upper abdominal pain location'],
            'Score': 0
        },
        'C. Splenic laceration': {
            'Consistent': ['Difficult colonoscopy (traction on splenocolic ligament)', 'LUQ / Epigastric pain', 'Left-shoulder discomfort (Kehr\'s sign)', 'Tachycardia and hypotension (hemorrhagic shock)', 'Massive Hemoglobin Drop', 'LUQ tenderness, guarding, peritoneal signs'],
            'Inconsistent': [],
            'Score': 0
        },
        'D. Postpolypectomy syndrome': {
            'Consistent': ['Abdominal pain after colonoscopy'],
            'Inconsistent': ['No polypectomy was performed (rules out diagnosis)', 'Primarily inflammatory, not associated with massive hemorrhage'],
            'Score': 0
        }
    }

    print("Evaluating patient findings against potential diagnoses...\n")

    # Step 3: Score each diagnosis based on the evidence.
    # Positive points for consistent findings, high points for classic signs, and negative for contradictions.
    
    # Diagnosis A: Colonic Perforation
    diagnoses_features['A. Colonic perforation']['Score'] += 2 # Consistent with difficult colonoscopy and peritoneal signs
    print("Scoring 'A. Colonic perforation':")
    print("  - Peritoneal signs are present. (+1 point)")
    print("  - It is a known complication of colonoscopy. (+1 point)")
    print("  - However, the pain location (LUQ + shoulder) and massive bleed are more indicative of another cause.")
    
    # Diagnosis B: Lower GI Bleeding
    diagnoses_features['B. Lower GI bleeding']['Score'] += 1 # Consistent with bleed after colonoscopy
    print("\nScoring 'B. Lower GI bleeding':")
    print("  - Patient is clearly bleeding. (+1 point)")
    print("  - However, the pain is in the upper abdomen, not lower, and the significant pain itself is atypical for some common causes of LGI bleeds.")

    # Diagnosis D: Postpolypectomy Syndrome
    diagnoses_features['D. Postpolypectomy syndrome']['Score'] -= 10 # Ruled out by a key fact.
    print("\nScoring 'D. Postpolypectomy syndrome':")
    print("  - The case explicitly states 'No polypectomy was performed'. This diagnosis is ruled out. (-10 points)")

    # Diagnosis C: Splenic Laceration
    # This matches all key high-yield findings.
    score_c = 0
    reasoning_c = []
    if 'Difficult colonoscopy' in diagnoses_features['C. Splenic laceration']['Consistent']:
        score_c += 1
        reasoning_c.append("  - Plausible mechanism: difficult colonoscopy can cause traction on the spleen. (+1 point)")
    if 'LUQ / Epigastric pain' in diagnoses_features['C. Splenic laceration']['Consistent']:
        score_c += 2
        reasoning_c.append("  - Matches the classic location of pain for splenic injury. (+2 points)")
    if 'Left-shoulder discomfort (Kehr\'s sign)' in diagnoses_features['C. Splenic laceration']['Consistent']:
        score_c += 5
        reasoning_c.append("  - Presence of Kehr's sign is a hallmark of diaphragmatic irritation from splenic bleed. (+5 points)")
    if 'Tachycardia and hypotension (hemorrhagic shock)' in diagnoses_features['C. Splenic laceration']['Consistent']:
        score_c += 3
        reasoning_c.append("  - Explains the rapid hemodynamic deterioration. (+3 points)")
    if 'Massive Hemoglobin Drop' in diagnoses_features['C. Splenic laceration']['Consistent']:
        score_c += 3
        reasoning_c.append("  - Consistent with bleeding from a highly vascular solid organ like the spleen. (+3 points)")

    diagnoses_features['C. Splenic laceration']['Score'] = score_c
    print("\nScoring 'C. Splenic laceration':")
    for r in reasoning_c:
        print(r)
    
    # Step 4: Display the final scores and the conclusion.
    print("\n--- Final Score Calculation ---")
    final_scores = {name: data['Score'] for name, data in diagnoses_features.items()}
    for name, score in final_scores.items():
        print(f"Diagnosis '{name}': Final Score = {score}")

    winner = max(final_scores, key=final_scores.get)
    print("\n--- Conclusion ---")
    print(f"The diagnosis with the highest score is '{winner}'.")
    print("\nThe combination of a difficult colonoscopy, acute onset of LUQ pain, referred left-shoulder pain (Kehr's sign), and rapid progression to hemorrhagic shock is the classic presentation for a splenic laceration.")

solve_medical_case()
<<<C>>>