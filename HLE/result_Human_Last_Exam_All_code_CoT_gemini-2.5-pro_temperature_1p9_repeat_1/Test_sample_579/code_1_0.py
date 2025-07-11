def solve_medical_case():
    """
    Analyzes a clinical case to determine the most likely diagnosis.
    This function simulates a diagnostic process by scoring potential conditions
    based on how well they match the patient's symptoms and risk factors.
    """

    # --- Patient Data from the Vignette ---
    patient_profile = {
        'risk_factors': {'obesity', 'smoking'},
        'locations': {'axillary', 'inframammary', 'inguinal'},
        'lesions': {'bullae', 'plaques', 'purulent_nodules'}
    }

    # --- Database of Diagnostic Criteria ---
    diagnoses = {
        'A. Malignant Intertrigo': {'score': 0, 'reasoning': []},
        'B. Allergic contact dermatitis': {'score': 0, 'reasoning': []},
        'C. Hidradenitis Suppurativa': {'score': 0, 'reasoning': []},
        'D. Atopic dermatitis': {'score': 0, 'reasoning': []},
        'E. Psoriasis': {'score': 0, 'reasoning': []}
    }

    # --- Scoring Logic ---

    # C. Hidradenitis Suppurativa (HS)
    # Strong association with obesity and smoking
    if 'obesity' in patient_profile['risk_factors']:
        diagnoses['C. Hidradenitis Suppurativa']['score'] += 2
        diagnoses['C. Hidradenitis Suppurativa']['reasoning'].append(('Obesity Risk Factor', 2))
    if 'smoking' in patient_profile['risk_factors']:
        diagnoses['C. Hidradenitis Suppurativa']['score'] += 2
        diagnoses['C. Hidradenitis Suppurativa']['reasoning'].append(('Smoking Risk Factor', 2))
    # Classic locations for HS
    location_match_hs = len(patient_profile['locations'].intersection({'axillary', 'inguinal', 'inframammary'}))
    diagnoses['C. Hidradenitis Suppurativa']['score'] += location_match_hs * 2
    diagnoses['C. Hidradenitis Suppurativa']['reasoning'].append(('Classic Locations', location_match_hs * 2))
    # Purulent nodules are a hallmark of HS
    if 'purulent_nodules' in patient_profile['lesions']:
        diagnoses['C. Hidradenitis Suppurativa']['score'] += 4
        diagnoses['C. Hidradenitis Suppurativa']['reasoning'].append(('Purulent Nodules (Hallmark)', 4))
    
    # E. Psoriasis (Inverse)
    if 'obesity' in patient_profile['risk_factors']:
        diagnoses['E. Psoriasis']['score'] += 1
        diagnoses['E. Psoriasis']['reasoning'].append(('Obesity Risk Factor', 1))
    # Plaques and locations fit inverse psoriasis
    if 'plaques' in patient_profile['lesions']:
        location_match_psoriasis = len(patient_profile['locations'].intersection({'axillary', 'inguinal', 'inframammary'}))
        diagnoses['E. Psoriasis']['score'] += 1 * location_match_psoriasis
        diagnoses['E. Psoriasis']['reasoning'].append(('Plaques in Folds', 1 * location_match_psoriasis))
    # Purulent nodules are a strong negative sign for psoriasis
    if 'purulent_nodules' in patient_profile['lesions']:
        diagnoses['E. Psoriasis']['score'] -= 3
        diagnoses['E. Psoriasis']['reasoning'].append(('Purulent Nodules (Negative Sign)', -3))

    # --- Find the best diagnosis ---
    best_diagnosis = max(diagnoses, key=lambda d: diagnoses[d]['score'])
    
    # --- Output the Reasoning ---
    print("Clinical Reasoning Analysis:")
    print("-" * 30)
    for diagnosis, data in sorted(diagnoses.items(), key=lambda item: item[1]['score'], reverse=True):
        print(f"Diagnosis: {diagnosis}\nFinal Score: {data['score']}\n")

    # --- Print the detailed equation for the top diagnosis ---
    top_dx_data = diagnoses[best_diagnosis]
    equation_parts = [f"{reason} ({score})" for reason, score in top_dx_data['reasoning']]
    equation_str = " + ".join(equation_parts).replace("+ -", "- ")
    
    print("=" * 30)
    print(f"Detailed Calculation for the Most Likely Diagnosis:\n")
    print(f"Diagnosis: {best_diagnosis}")
    # This line prints the numbers in the final equation as requested.
    print(f"Score Calculation: {equation_str} = {top_dx_data['score']}")
    print("=" * 30)
    
    print("\nConclusion:")
    print("The patient presents with obesity and a smoking history, which are major risk factors for Hidradenitis Suppurativa (HS). The involvement of the axillary, inframammary, and inguinal folds is classic for HS. Most importantly, the presence of purulent nodules is a hallmark feature of HS, making it the most likely diagnosis over others like Psoriasis, which does not cause deep, purulent nodules.")

solve_medical_case()
<<<C>>>