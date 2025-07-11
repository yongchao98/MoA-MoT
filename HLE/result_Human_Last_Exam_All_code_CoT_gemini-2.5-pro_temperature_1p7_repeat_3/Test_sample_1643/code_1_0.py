import sys

def diagnose_antiphospholipid_syndrome():
    """
    Evaluates if a patient meets the diagnostic criteria for Antiphospholipid Syndrome (APS)
    based on the revised Sapporo (Sydney) classification criteria.
    """

    # --- Patient Data ---
    # Clinical Data
    vte_events = 3

    # Lab Data from 3 months ago (T1)
    t1_antib2gp1_igm = 41  # N < 20 UI/L
    t1_ac_igm = 32         # N < 20 UI/L
    t1_ptt_la_ratio = 1.19 # N < 1.18
    t1_drvvt_ratio = 1.44  # N < 1.2

    # Lab Data from today (T2, >12 weeks later)
    t2_antib2gp1_igg = 21  # N < 20 UI/L
    t2_antib2gp1_igm = 29  # N < 20 UI/L
    t2_ac_igm = 47         # N < 20 UI/L
    t2_ptt_la_ratio = 1.17 # N < 1.18
    t2_drvvt_ratio = 1.51  # N < 1.2

    # --- Criteria Evaluation ---
    print("Evaluating patient for Antiphospholipid Syndrome (APS)")
    print("="*50)

    # 1. Clinical Criteria Assessment
    print("Step 1: Checking Clinical Criteria...")
    clinical_criterion_met = False
    if vte_events >= 1:
        print(f"- Patient has a history of {vte_events} venous thromboembolic (VTE) events.")
        print("- RESULT: The clinical criterion (Vascular Thrombosis) is MET.")
        clinical_criterion_met = True
    else:
        print("- Patient has no history of VTE.")
        print("- RESULT: The clinical criterion is NOT MET.")
    print("-" * 50)

    # 2. Laboratory Criteria Assessment
    print("Step 2: Checking Laboratory Criteria (must be positive on 2 occasions >12 weeks apart)...")
    lab_criterion_met = False
    
    # Lupus Anticoagulant (LA)
    la_positive = False
    # dRVVT is the more specific test and is strongly positive twice.
    if t1_drvvt_ratio > 1.2 and t2_drvvt_ratio > 1.2:
        print(f"- Lupus Anticoagulant (LA) assessment:")
        print(f"  - dRVVT ratio was {t1_drvvt_ratio} (T1) and {t2_drvvt_ratio} (T2). Both are > 1.2.")
        print("  - The PTT-LA ratio was {t1_ptt_la_ratio} (T1, positive) and {t2_ptt_la_ratio} (T2, negative).")
        print(f"  - RESULT: LA is considered persistently positive based on dRVVT. Criterion MET.")
        la_positive = True
        lab_criterion_met = True
    else:
        print("- RESULT: Lupus Anticoagulant is NOT persistently positive.")

    # Anticardiolipin (aCL) antibodies
    acl_positive = False
    if (t1_ac_igm > 20 and t2_ac_igm > 20):
         print(f"- Anticardiolipin (aCL) antibody assessment:")
         print(f"  - aCL IgM was {t1_ac_igm} UI/L (T1) and {t2_ac_igm} UI/L (T2). Both are > 20 UI/L (positive).")
         # Check for medium/high titer (>40)
         if t1_ac_igm > 40 or t2_ac_igm > 40:
             print(f"  - Titer reached medium/high level at T2 ({t2_ac_igm} UI/L).")
         print("  - RESULT: aCL IgM is persistently positive. Criterion MET.")
         acl_positive = True
         lab_criterion_met = True
    else:
        print("- RESULT: Anticardiolipin antibodies are NOT persistently positive.")

    # Anti-ß2-glycoprotein-I (anti-ß2GPI) antibodies
    b2gpi_positive = False
    if (t1_antib2gp1_igm > 20) and (t2_antib2gp1_igm > 20 or t2_antib2gp1_igg > 20):
        print(f"- Anti-ß2GP1 antibody assessment:")
        print(f"  - Anti-ß2GP1 IgM was {t1_antib2gp1_igm} UI/L at T1 (positive).")
        print(f"  - At T2, Anti-ß2GP1 IgM was {t2_antib2gp1_igm} UI/L (positive) and IgG was {t2_antib2gp1_igg} UI/L (positive).")
        print("  - RESULT: Anti-ß2GP1 antibodies are persistently positive. Criterion MET.")
        b2gpi_positive = True
        lab_criterion_met = True
    else:
        print("- RESULT: Anti-ß2GP1 antibodies are NOT persistently positive.")
        
    print("-" * 50)
    
    # 3. Final Diagnosis
    print("Step 3: Final Conclusion...")
    final_answer = ""
    if clinical_criterion_met and lab_criterion_met:
        print("The patient meets at least one clinical criterion AND at least one laboratory criterion.")
        final_answer = "Yes"
    else:
        print("The patient does not meet both the clinical and laboratory criteria for APS.")
        final_answer = "No"

    print(f"\nDoes this patient categorize as having antiphospholipid syndrome? {final_answer}")
    
    # Print final answer in specified format
    # Redirecting to stdout to ensure it is the last line.
    sys.stdout.write(f'<<<{final_answer}>>>\n')

diagnose_antiphospholipid_syndrome()