import sys

def diagnose_aps():
    """
    Analyzes patient data against APS criteria and provides a diagnosis.
    """
    # --- Patient Data ---
    patient_data = {
        'vte_events': 3,
        'medications': ['DESOGESTREL 75µg', 'RIVAROXABAN 20mg'],
        'lab_3_months_ago': {
            'anti_b2gp1_igm': 41, 'anti_b2gp1_igg': 18,
            'acl_igm': 32, 'acl_igg': 9,
            'ptt_la_ratio': 1.19, 'drvv_ratio': 1.44
        },
        'lab_today': {
            'anti_b2gp1_igm': 29, 'anti_b2gp1_igg': 21,
            'acl_igm': 47, 'acl_igg': 7,
            'ptt_la_ratio': 1.17, 'drvv_ratio': 1.51
        },
        'normals': {
            'antibody': 20, # for aCL and anti-b2gp1
            'ptt_la_ratio': 1.18,
            'drvv_ratio': 1.2
        }
    }

    print("--- Antiphospholipid Syndrome (APS) Diagnostic Check ---")
    
    # --- 1. Check Clinical Criteria ---
    print("\nStep 1: Evaluating Clinical Criteria...")
    vte_events = patient_data['vte_events']
    clinical_criterion_met = vte_events >= 1
    print(f"Does patient meet vascular thrombosis criterion (>= 1 event)?")
    print(f"Calculation: Patient has {vte_events} VTE events, which is >= 1.")
    print(f"Result: {'Yes' if clinical_criterion_met else 'No'}")

    # --- 2. Check Laboratory Criteria ---
    print("\nStep 2: Evaluating Laboratory Criteria (Persistent Positivity >12 weeks apart)...")
    
    lab1 = patient_data['lab_3_months_ago']
    lab2 = patient_data['lab_today']
    normals = patient_data['normals']
    
    lab_markers_positive = []

    # Check aCL IgM
    acl_igm_p = lab1['acl_igm'] > normals['antibody'] and lab2['acl_igm'] > normals['antibody']
    if acl_igm_p: lab_markers_positive.append("aCL IgM")
    print(f"Is aCL IgM persistently positive? (Test 1: {lab1['acl_igm']}, Test 2: {lab2['acl_igm']}, Normal < {normals['antibody']}) -> {'Yes' if acl_igm_p else 'No'}")

    # Check anti-b2gp1 IgM
    b2gp1_igm_p = lab1['anti_b2gp1_igm'] > normals['antibody'] and lab2['anti_b2gp1_igm'] > normals['antibody']
    if b2gp1_igm_p: lab_markers_positive.append("anti-b2gp1 IgM")
    print(f"Is anti-b2gp1 IgM persistently positive? (Test 1: {lab1['anti_b2gp1_igm']}, Test 2: {lab2['anti_b2gp1_igm']}, Normal < {normals['antibody']}) -> {'Yes' if b2gp1_igm_p else 'No'}")

    # Check dRVVT
    drvv_p = lab1['drvv_ratio'] > normals['drvv_ratio'] and lab2['drvv_ratio'] > normals['drvv_ratio']
    if drvv_p: lab_markers_positive.append("dRVVT (Lupus Anticoagulant)")
    print(f"Is dRVVT persistently positive? (Test 1: {lab1['drvv_ratio']}, Test 2: {lab2['drvv_ratio']}, Normal < {normals['drvv_ratio']}) -> {'Yes' if drvv_p else 'No'}")
    if drvv_p and 'RIVAROXABAN 20mg' in patient_data['medications']:
        print("  -> WARNING: Positive dRVVT is likely a false positive due to Rivaroxaban.")

    laboratory_criterion_met = len(lab_markers_positive) > 0
    print(f"\nDoes patient meet laboratory criterion (>= 1 persistent marker)?")
    if laboratory_criterion_met:
        print(f"Result: Yes. Persistently positive markers found: {', '.join(lab_markers_positive)}")
    else:
        print("Result: No. No persistently positive markers found.")

    # --- 3. Final Conclusion ---
    print("\n--- Final Diagnosis ---")
    if clinical_criterion_met and laboratory_criterion_met:
        final_answer = "Yes"
        print("Conclusion: Patient meets BOTH clinical and laboratory criteria.")
    else:
        final_answer = "No"
        print("Conclusion: Patient does not meet both required criteria for an APS diagnosis.")
        if not clinical_criterion_met: print("- Clinical criterion NOT met.")
        if not laboratory_criterion_met: print("- Laboratory criterion NOT met.")
    
    # Suppress final answer from print to stdout to match specified format
    original_stdout = sys.stdout 
    sys.stdout = open('/dev/null', 'w')
    print(f"<<<{final_answer}>>>")
    sys.stdout.close()
    sys.stdout = original_stdout
    
if __name__ == '__main__':
    diagnose_aps()
