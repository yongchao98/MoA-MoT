import sys

# Redirect print to a string to capture it for the final output format if needed,
# though for this task, direct printing is fine.

def check_aps_criteria():
    """
    Evaluates a patient's data against the APS Sydney criteria.
    """
    # --- Patient Data ---
    # Clinical
    vte_events = 3
    
    # Lab Results (T1 = 3 months ago, T2 = today)
    lab_t1 = {
        "aCL_IgM": 41, "aCL_IgG": 9,
        "ab2GPI_IgM": 32, "ab2GPI_IgG": 18,
        "dRVVT_ratio": 1.44
    }
    lab_t2 = {
        "aCL_IgM": 47, "aCL_IgG": 7,
        "ab2GPI_IgM": 29, "ab2GPI_IgG": 21,
        "dRVVT_ratio": 1.51
    }
    
    # --- Criteria Cutoffs (based on provided normals) ---
    aCL_cutoff = 20
    ab2GPI_cutoff = 20
    dRVVT_cutoff = 1.2
    # Medium/high titer for aCL is often defined as >40 GPL/MPL
    aCL_medium_titer_cutoff = 40

    print("Evaluating patient for Antiphospholipid Syndrome (APS) based on Sydney Criteria\n")

    # 1. Evaluate Clinical Criteria
    print("--- Step 1: Clinical Criteria Evaluation ---")
    clinical_criterion_met = False
    if vte_events >= 1:
        clinical_criterion_met = True
        print(f"Criterion: At least 1 episode of venous or arterial thrombosis.")
        print(f"Result: Patient has had {vte_events} VTE events. Clinical criterion is MET.\n")
    else:
        # Also check for pregnancy morbidity, but VTE is sufficient
        print("Result: No definite history of vascular thrombosis. Clinical criterion is NOT MET.\n")

    # 2. Evaluate Laboratory Criteria (must be persistent over >= 12 weeks)
    print("--- Step 2: Laboratory Criteria Evaluation ---")
    
    # Lupus Anticoagulant (LA) via dRVVT
    la_positive_t1 = lab_t1['dRVVT_ratio'] > dRVVT_cutoff
    la_positive_t2 = lab_t2['dRVVT_ratio'] > dRVVT_cutoff
    la_criterion_met = la_positive_t1 and la_positive_t2
    print("Criterion A: Lupus Anticoagulant (LA) positive on 2 occasions >= 12 weeks apart.")
    print(f" - Test 1 (3 months ago): dRVVT ratio was {lab_t1['dRVVT_ratio']} (Normal < {dRVVT_cutoff}). Positive: {la_positive_t1}")
    print(f" - Test 2 (Today): dRVVT ratio is {lab_t2['dRVVT_ratio']} (Normal < {dRVVT_cutoff}). Positive: {la_positive_t2}")
    print(f"Result: LA criterion is {'MET' if la_criterion_met else 'NOT MET'}.\n")

    # Anticardiolipin (aCL) antibody
    # Must be medium-high titer (>40) or persistently positive
    acl_positive_t1 = lab_t1['aCL_IgM'] > aCL_cutoff or lab_t1['aCL_IgG'] > aCL_cutoff
    acl_positive_t2 = lab_t2['aCL_IgM'] > aCL_cutoff or lab_t2['aCL_IgG'] > aCL_cutoff
    acl_medium_titer_present = lab_t1['aCL_IgM'] > aCL_medium_titer_cutoff or lab_t2['aCL_IgM'] > aCL_medium_titer_cutoff
    acl_criterion_met = (acl_positive_t1 and acl_positive_t2) and acl_medium_titer_present
    print("Criterion B: Anticardiolipin (aCL) Ab (IgG/IgM) positive on 2 occasions >= 12 weeks apart.")
    print(f" - Test 1 (3 months ago): aCL IgM was {lab_t1['aCL_IgM']} (Normal < {aCL_cutoff}). Positive: {acl_positive_t1}")
    print(f" - Test 2 (Today): aCL IgM is {lab_t2['aCL_IgM']} (Normal < {aCL_cutoff}). Positive: {acl_positive_t2}")
    print(f" - Additional check: At least one value is medium-titer (> {aCL_medium_titer_cutoff}). Status: {acl_medium_titer_present}")
    print(f"Result: aCL criterion is {'MET' if acl_criterion_met else 'NOT MET'}.\n")
    
    # Anti-β2 glycoprotein I (anti-β2GPI) antibody
    ab2gpi_positive_t1 = lab_t1['ab2GPI_IgM'] > ab2GPI_cutoff or lab_t1['ab2GPI_IgG'] > ab2GPI_cutoff
    ab2gpi_positive_t2 = lab_t2['ab2GPI_IgM'] > ab2GPI_cutoff or lab_t2['ab2GPI_IgG'] > ab2GPI_cutoff
    ab2gpi_criterion_met = ab2gpi_positive_t1 and ab2gpi_positive_t2
    print("Criterion C: Anti-β2GPI Ab (IgG/IgM) positive on 2 occasions >= 12 weeks apart.")
    print(f" - Test 1 (3 months ago): anti-β2GPI IgM was {lab_t1['ab2GPI_IgM']} (Normal < {ab2GPI_cutoff}). Positive: {ab2gpi_positive_t1}")
    print(f" - Test 2 (Today): anti-β2GPI IgM is {lab_t2['ab2GPI_IgM']} and IgG is {lab_t2['ab2GPI_IgG']} (Normal < {ab2GPI_cutoff}). Positive: {ab2gpi_positive_t2}")
    print(f"Result: anti-β2GPI criterion is {'MET' if ab2gpi_criterion_met else 'NOT MET'}.\n")

    laboratory_criterion_met = la_criterion_met or acl_criterion_met or ab2gpi_criterion_met
    print(f"Overall Laboratory Criteria Status: {'MET' if laboratory_criterion_met else 'NOT MET'} (at least one of A, B, or C must be met).\n")

    # 3. Final Conclusion
    print("--- Step 3: Final Diagnosis ---")
    final_diagnosis_met = clinical_criterion_met and laboratory_criterion_met
    final_answer = "Yes" if final_diagnosis_met else "No"
    
    print(f"Diagnosis of APS requires: (Clinical Criterion MET) AND (Laboratory Criterion MET)")
    print(f"Patient's status: ({clinical_criterion_met}) AND ({laboratory_criterion_met})")
    print(f"Conclusion: The patient categorizes as having antiphospholipid syndrome.")
    
    return final_answer

if __name__ == '__main__':
    answer = check_aps_criteria()
    print(f"<<<{answer}>>>")