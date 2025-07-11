import sys

# Disable writing of bytecode files
sys.dont_write_bytecode = True

def diagnose_aps():
    """
    Analyzes patient data against the Sydney criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    # Clinical
    vte_events = 3
    # Lab 1 (3 months ago)
    lab1_dRVVT_ratio = 1.44
    lab1_aCL_IgM = 32
    lab1_ab2GP1_IgM = 41
    lab1_ab2GP1_IgG = 18
    # Lab 2 (Today)
    lab2_dRVVT_ratio = 1.51
    lab2_aCL_IgM = 47
    lab2_ab2GP1_IgM = 29
    lab2_ab2GP1_IgG = 21

    # Normal Values
    dRVVT_norm = 1.2
    aCL_ab2GP1_norm = 20

    clinical_criteria_met = False
    lab_criteria_met = False

    print("Step 1: Assessing Clinical Criteria for APS")
    print("-------------------------------------------")
    print(f"The patient has a history of {vte_events} venous thromboembolism (VTE) events.")
    if vte_events >= 1:
        clinical_criteria_met = True
        print("Conclusion: At least one documented vascular thrombosis event has occurred. The clinical criterion is MET.")
    else:
        print("Conclusion: The clinical criterion is NOT MET.")
    print("\n")


    print("Step 2: Assessing Laboratory Criteria for APS (Persistence over >12 weeks)")
    print("-------------------------------------------------------------------")
    
    # Check for Lupus Anticoagulant (LA)
    la_positive_1 = lab1_dRVVT_ratio > dRVVT_norm
    la_positive_2 = lab2_dRVVT_ratio > dRVVT_norm
    persistent_la = la_positive_1 and la_positive_2
    print(f"Lupus Anticoagulant (dRVVT):")
    print(f" - 3 months ago: dRVVT ratio was {lab1_dRVVT_ratio} (Normal < {dRVVT_norm}) -> {'Positive' if la_positive_1 else 'Negative'}")
    print(f" - Today: dRVVT ratio was {lab2_dRVVT_ratio} (Normal < {dRVVT_norm}) -> {'Positive' if la_positive_2 else 'Negative'}")
    if persistent_la:
        print(" -> Finding: Persistent Lupus Anticoagulant detected.")
        print(" -> NOTE: Rivaroxaban (a DOAC) can cause false-positive dRVVT results. However, this is a positive finding based on the provided numbers.")
        lab_criteria_met = True
    else:
        print(" -> Finding: Lupus Anticoagulant is not persistently positive based on dRVVT.")
    print("")

    # Check for anticardiolipin (aCL) antibodies
    acl_positive_1 = lab1_aCL_IgM > aCL_ab2GP1_norm
    acl_positive_2 = lab2_aCL_IgM > aCL_ab2GP1_norm
    persistent_acl = acl_positive_1 and acl_positive_2
    print("Anticardiolipin Antibodies (aCL IgM):")
    print(f" - 3 months ago: aCL IgM was {lab1_aCL_IgM} UI/L (Normal < {aCL_ab2GP1_norm}) -> {'Positive' if acl_positive_1 else 'Negative'}")
    print(f" - Today: aCL IgM was {lab2_aCL_IgM} UI/L (Normal < {aCL_ab2GP1_norm}) -> {'Positive' if acl_positive_2 else 'Negative'}")
    if persistent_acl:
        print(" -> Finding: Persistent aCL IgM antibodies detected.")
        lab_criteria_met = True
    else:
        print(" -> Finding: aCL IgM antibodies are not persistently positive.")
    print("")

    # Check for anti-B2GP1 antibodies
    ab2gp1_positive_1 = lab1_ab2GP1_IgM > aCL_ab2GP1_norm or lab1_ab2GP1_IgG > aCL_ab2GP1_norm
    ab2gp1_positive_2 = lab2_ab2GP1_IgM > aCL_ab2GP1_norm or lab2_ab2GP1_IgG > aCL_ab2GP1_norm
    persistent_ab2gp1 = ab2gp1_positive_1 and ab2gp1_positive_2
    print("Anti-ß2GP1 Antibodies:")
    print(f" - 3 months ago: IgM was {lab1_ab2GP1_IgM}, IgG was {lab1_ab2GP1_IgG} (Normal < {aCL_ab2GP1_norm}) -> {'Positive' if ab2gp1_positive_1 else 'Negative'}")
    print(f" - Today: IgM was {lab2_ab2GP1_IgM}, IgG was {lab2_ab2GP1_IgG} (Normal < {aCL_ab2GP1_norm}) -> {'Positive' if ab2gp1_positive_2 else 'Negative'}")
    if persistent_ab2gp1:
        print(" -> Finding: Persistent anti-ß2GP1 antibodies detected.")
        lab_criteria_met = True
    else:
        print(" -> Finding: Anti-ß2GP1 antibodies are not persistently positive.")
    print("")

    if lab_criteria_met:
        print("Conclusion: At least one laboratory criterion (persistently positive antiphospholipid antibodies) is MET.")
    else:
        print("Conclusion: The laboratory criteria are NOT MET.")
    print("\n")


    print("Step 3: Final Diagnosis")
    print("------------------------")
    final_answer = ""
    if clinical_criteria_met and lab_criteria_met:
        final_answer = "Yes"
        print("The patient meets at least one clinical criterion (Vascular Thrombosis) AND at least one laboratory criterion (Persistent antiphospholipid antibodies).")
        print("Therefore, the patient categorizes as having antiphospholipid syndrome.")
    else:
        final_answer = "No"
        print("The patient does not meet both the clinical and laboratory criteria required for a diagnosis of antiphospholipid syndrome.")
    
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    diagnose_aps()