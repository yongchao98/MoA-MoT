def diagnose_case():
    """
    This function analyzes clinical findings from a case study to determine the most likely diagnosis.
    It uses a scoring system to weigh the evidence for each possibility.
    """

    # --- Clinical Findings from the Case ---
    # Finding: (Description, Weight for its significance)
    findings = {
        'procedure_difficulty': ('Difficult Colonoscopy', 1),
        'no_polypectomy': ('No Polypectomy Performed', 5), # Very high weight as it can rule out a diagnosis
        'pain_location': ('Left Upper Quadrant Pain', 3), # Highly specific
        'referred_pain': ("Left Shoulder Pain (Kehr's Sign)", 4), # Highly specific
        'hemorrhage_signs': ('Rapid Massive Hemoglobin Drop & Shock', 4), # Hallmark of the problem
        'peritoneal_signs': ('Abdominal Distension and Guarding', 2)
    }

    # --- Potential Diagnoses ---
    diagnoses = {
        'A. Colonic perforation': 0,
        'B. Lower GI bleeding': 0,
        'C. Splenic laceration': 0,
        'D. Postpolypectomy syndrome': 0
    }

    # --- Scoring Logic ---

    # Score for 'C. Splenic laceration'
    score_c_calc = []
    score_c = 0
    # Strong association with all major findings
    score_c += findings['procedure_difficulty'][1]; score_c_calc.append(str(findings['procedure_difficulty'][1]))
    score_c += findings['pain_location'][1]; score_c_calc.append(str(findings['pain_location'][1]))
    score_c += findings['referred_pain'][1]; score_c_calc.append(str(findings['referred_pain'][1]))
    score_c += findings['hemorrhage_signs'][1]; score_c_calc.append(str(findings['hemorrhage_signs'][1]))
    score_c += findings['peritoneal_signs'][1]; score_c_calc.append(str(findings['peritoneal_signs'][1]))
    diagnoses['C. Splenic laceration'] = score_c

    # Score for 'A. Colonic perforation'
    score_a_calc = []
    score_a = 0
    # Possible but less likely fit
    score_a += findings['procedure_difficulty'][1]; score_a_calc.append(str(findings['procedure_difficulty'][1]))
    score_a += findings['peritoneal_signs'][1]; score_a_calc.append(str(findings['peritoneal_signs'][1]))
    # Hemorrhage of this magnitude is less typical than for splenic injury
    score_a += 1; score_a_calc.append("1") # Points for some bleeding, but not full weight
    diagnoses['A. Colonic perforation'] = score_a

    # Score for 'B. Lower GI bleeding'
    score_b_calc = []
    score_b = 0
    # Unlikely due to pain location and no polypectomy
    score_b += findings['hemorrhage_signs'][1] / 2; score_b_calc.append(str(findings['hemorrhage_signs'][1]/2)) # Bleeding is present, but signs point away from this source
    score_b -= findings['no_polypectomy'][1]; score_b_calc.append(str(-findings['no_polypectomy'][1])) # Negative points
    diagnoses['B. Lower GI bleeding'] = score_b

    # Score for 'D. Postpolypectomy syndrome'
    score_d_calc = []
    score_d = 0
    # Ruled out by no polypectomy
    score_d -= findings['no_polypectomy'][1]; score_d_calc.append(str(-findings['no_polypectomy'][1])) # Negative points
    diagnoses['D. Postpolypectomy syndrome'] = score_d

    # --- Print Results ---
    print("Clinical Reasoning Score Sheet:\n")
    print(f"A. Colonic perforation score calculation: {' + '.join(score_a_calc)} = {diagnoses['A. Colonic perforation']}")
    print(f"B. Lower GI bleeding score calculation: {' + '.join(score_b_calc)} = {diagnoses['B. Lower GI bleeding']}")
    print(f"C. Splenic laceration score calculation: {' + '.join(score_c_calc)} = {diagnoses['C. Splenic laceration']}")
    print(f"D. Postpolypectomy syndrome score calculation: {' + '.join(score_d_calc)} = {diagnoses['D. Postpolypectomy syndrome']}")
    print("-" * 30)

    # Determine the most likely diagnosis
    most_likely_diagnosis = max(diagnoses, key=diagnoses.get)
    print(f"The most likely diagnosis is: {most_likely_diagnosis}")
    print("\nThis diagnosis best explains the combination of a difficult colonoscopy, left upper quadrant pain, left shoulder pain (Kehr's sign), and massive intra-abdominal hemorrhage.")

# Execute the diagnostic analysis
diagnose_case()