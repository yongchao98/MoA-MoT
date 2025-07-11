def diagnose_antiphospholipid_syndrome():
    """
    Analyzes patient data against the APS diagnostic criteria.
    The international criteria require at least ONE clinical criterion AND ONE laboratory criterion.
    The lab criterion must be persistent, i.e., present on two or more occasions at least 12 weeks apart.
    """

    # --- Patient Data ---
    patient_history = {
        "vte_events": 3,
        "pregnancy_with_pe": True
    }

    labs_3_months_ago = {
        "anti_b2gp1_igm": 41, "anti_b2gp1_igg": 18,
        "anticardiolipin_igm": 32, "anticardiolipin_igg": 9,
        "ptt_la_ratio": 1.19, "drvv_ratio": 1.44
    }

    labs_today = {
        "anti_b2gp1_igm": 29, "anti_b2gp1_igg": 21,
        "anticardiolipin_igm": 47, "anticardiolipin_igg": 7,
        "ptt_la_ratio": 1.17, "drvv_ratio": 1.51
    }

    # --- Criteria Normal Values ---
    normals = {
        "antibody_level": 20,
        "ptt_la_ratio": 1.18,
        "drvv_ratio": 1.2
    }
    
    # --- Analysis ---
    clinical_criterion_met = False
    lab_criterion_met = False
    
    # Step 1: Check Clinical Criteria
    print("Step 1: Evaluating Clinical Criteria...")
    if patient_history["vte_events"] >= 1:
        clinical_criterion_met = True
        print(f"-> Clinical Criterion MET: Patient has a history of {patient_history['vte_events']} vascular thrombotic events.")
    else:
        print("-> Clinical Criterion NOT MET.")
    print("-" * 50)

    # Step 2: Check Laboratory Criteria (Persistence over >= 12 weeks)
    print("Step 2: Evaluating Laboratory Criteria...")
    
    # Check for persistent Lupus Anticoagulant (LA)
    la_positive_t1 = labs_3_months_ago["drvv_ratio"] > normals["drvv_ratio"] or labs_3_months_ago["ptt_la_ratio"] > normals["ptt_la_ratio"]
    la_positive_t2 = labs_today["drvv_ratio"] > normals["drvv_ratio"] or labs_today["ptt_la_ratio"] > normals["ptt_la_ratio"]
    if la_positive_t1 and la_positive_t2:
        lab_criterion_met = True
        print(f"-> Lab Criterion MET: Persistent Lupus Anticoagulant detected.")
        print(f"   - 3 months ago: dRVVT ratio was {labs_3_months_ago['drvv_ratio']} (Normal < {normals['drvv_ratio']})")
        print(f"   - Today: dRVVT ratio is {labs_today['drvv_ratio']} (Normal < {normals['drvv_ratio']})")

    # Check for persistent Anticardiolipin (aCL) antibodies
    acl_igm_persistent = labs_3_months_ago["anticardiolipin_igm"] > normals["antibody_level"] and labs_today["anticardiolipin_igm"] > normals["antibody_level"]
    acl_igg_persistent = labs_3_months_ago["anticardiolipin_igg"] > normals["antibody_level"] and labs_today["anticardiolipin_igg"] > normals["antibody_level"]
    if acl_igm_persistent or acl_igg_persistent:
        lab_criterion_met = True
        print(f"-> Lab Criterion MET: Persistent Anticardiolipin antibodies detected.")
        print(f"   - aCL IgM was {labs_3_months_ago['anticardiolipin_igm']} and is now {labs_today['anticardiolipin_igm']} (Normal < {normals['antibody_level']})")

    # Check for persistent Anti-B2GP1 antibodies
    ab2_igm_persistent = labs_3_months_ago["anti_b2gp1_igm"] > normals["antibody_level"] and labs_today["anti_b2gp1_igm"] > normals["antibody_level"]
    ab2_igg_persistent = labs_3_months_ago["anti_b2gp1_igg"] > normals["antibody_level"] and labs_today["anti_b2gp1_igg"] > normals["antibody_level"]
    if ab2_igm_persistent or ab2_igg_persistent:
        lab_criterion_met = True
        print(f"-> Lab Criterion MET: Persistent Anti-beta-2-glycoprotein-I antibodies detected.")
        print(f"   - anti-B2GP1 IgM was {labs_3_months_ago['anti_b2gp1_igm']} and is now {labs_today['anti_b2gp1_igm']} (Normal < {normals['antibody_level']})")
        
    if not lab_criterion_met:
        print("-> No persistent laboratory marker found. Lab Criterion NOT MET.")
    print("-" * 50)

    # Step 3: Final Diagnosis
    print("Step 3: Final Conclusion...")
    if clinical_criterion_met and lab_criterion_met:
        final_answer = "Yes"
        print("-> RESULT: The patient meets BOTH clinical and laboratory criteria.")
    else:
        final_answer = "No"
        print("-> RESULT: The patient does NOT meet the full criteria for Antiphospholipid Syndrome.")

    return final_answer

if __name__ == "__main__":
    answer = diagnose_antiphospholipid_syndrome()
    print(f"<<<{answer}>>>")