def check_nmr_signal_count():
    """
    This function verifies the number of distinct 1H NMR signals for the final product,
    Ethyl 1-cyanocyclohexanecarboxylate.

    The logic is broken down into analyzing the two main parts of the molecule:
    1. The ethyl ester group (-COOCH2CH3).
    2. The substituted cyclohexane ring.
    """

    # --- Part 1: Analysis of the Ethyl Group ---
    # The ethyl group consists of a methylene (-CH2-) and a methyl (-CH3-) group.
    # The three protons of the methyl group are chemically equivalent.
    # The two protons of the methylene group are also chemically equivalent.
    # These two groups are distinct from each other.
    signals_from_ethyl_group = 2  # One for -CH2-, one for -CH3-

    # --- Part 2: Analysis of the Cyclohexane Ring ---
    # The final product, Ethyl 1-cyanocyclohexanecarboxylate, has a plane of
    # symmetry that passes through C1 (the substituted carbon), C4, the cyano group,
    # and the ester group.
    # This symmetry makes pairs of carbons equivalent: C2 is equivalent to C6, and C3 is equivalent to C5.
    #
    # Crucially, C1 is a prochiral center because it is bonded to four different groups
    # (CN, COOEt, the C2/C6 path of the ring, and the C3/C5 path of the ring).
    # The presence of this prochiral center makes the two geminal protons on every
    # CH2 group in the ring diastereotopic (i.e., chemically non-equivalent).

    # Signals from C2 and C6:
    # C2 and C6 are equivalent by symmetry. However, on each of these carbons, the two
    # protons are diastereotopic. This results in 2 distinct signals for the 4 protons on C2/C6.
    signals_from_c2_c6 = 2

    # Signals from C3 and C5:
    # Similar to C2/C6, these positions are equivalent by symmetry, but the geminal
    # protons are diastereotopic. This results in 2 distinct signals for the 4 protons on C3/C5.
    signals_from_c3_c5 = 2

    # Signals from C4:
    # C4 is unique and lies on the plane of symmetry. Its two protons are diastereotopic
    # due to their relationship to the prochiral C1. This results in 2 distinct signals.
    signals_from_c4 = 2

    # Total signals from the ring. C1 has no protons.
    signals_from_ring = signals_from_c2_c6 + signals_from_c3_c5 + signals_from_c4

    # --- Part 3: Total Signal Count ---
    total_signals = signals_from_ethyl_group + signals_from_ring

    # --- Part 4: Verification ---
    # The question's option C corresponds to 8 signals.
    # The provided answer also concludes there are 8 signals.
    expected_answer = 8

    if total_signals == expected_answer:
        # The logical breakdown (2 signals from ethyl + 6 from ring) matches the expected total.
        return "Correct"
    else:
        reason = (
            f"The calculated total number of signals is {total_signals}, which does not match the expected answer of {expected_answer}.\n"
            f"Calculation breakdown:\n"
            f"- Signals from ethyl group: {signals_from_ethyl_group}\n"
            f"- Signals from cyclohexane ring: {signals_from_ring}\n"
            f"The reasoning in the provided answer leads to 8 signals, but the final conclusion might be wrong or there is a flaw in this verification code's logic."
        )
        return reason

# Execute the check
result = check_nmr_signal_count()
print(result)