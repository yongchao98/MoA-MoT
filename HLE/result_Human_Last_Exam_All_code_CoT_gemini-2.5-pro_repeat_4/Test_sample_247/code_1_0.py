def check_tcr_sequencing_reach(read_length, tcr_chain='beta'):
    """
    Calculates if a given sequencing read length can reach the CDR3 region
    from the 3' poly(A) tail of a TCR transcript.

    Args:
        read_length (int): The length of the sequencing read starting from the 3' end.
        tcr_chain (str): 'alpha' or 'beta' to specify the TCR chain.
    """

    # Approximate lengths of human TCR regions in base pairs (bp).
    # These are conservative estimates.
    # The 3' UTR (Untranslated Region) is highly variable in length.
    # The Constant region is downstream of the V(D)J (CDR3) region.
    UTR_3_PRIME_LENGTH_MIN = 150
    if tcr_chain == 'beta':
        # Human TRBC1/TRBC2 coding sequence length
        CONSTANT_REGION_LENGTH = 450
        chain_name = "TCR Beta"
    elif tcr_chain == 'alpha':
        # Human TRAC coding sequence length
        CONSTANT_REGION_LENGTH = 420
        chain_name = "TCR Alpha"
    else:
        print("Error: Invalid TCR chain specified. Use 'alpha' or 'beta'.")
        return

    # The distance from the start of the poly(A) tail to the CDR3 region
    # is at least the length of the 3' UTR plus the length of the constant region.
    distance_to_cdr3_min = UTR_3_PRIME_LENGTH_MIN + CONSTANT_REGION_LENGTH

    print(f"Analysis for {chain_name} Chain:")
    print(f"----------------------------------")
    print(f"Assumed 3' UTR length (minimum): {UTR_3_PRIME_LENGTH_MIN} bp")
    print(f"Assumed Constant Region length: {CONSTANT_REGION_LENGTH} bp")
    print(f"Sequencing Read 2 length: {read_length} bp")
    print(f"----------------------------------")
    print(f"Total distance from Poly(A) tail to start of CDR3 region (minimum estimate):")
    print(f"{UTR_3_PRIME_LENGTH_MIN} (3' UTR) + {CONSTANT_REGION_LENGTH} (Constant Region) = {distance_to_cdr3_min} bp")
    print(f"----------------------------------")

    if read_length >= distance_to_cdr3_min:
        print(f"Result: SUCCESS. The read length of {read_length} bp is likely sufficient to reach the CDR3 region.")
    else:
        print(f"Result: FAILURE. The read length of {read_length} bp is too short.")
        print(f"The read will end approximately {distance_to_cdr3_min - read_length} bp before reaching the CDR3 region.")
        print("\nConclusion: Priming from the 3' Poly(A) tail is not effective.")
        print("A different strategy, like using a primer in the Constant Region, is required.")


# The PhD student is using a 225 bp read (Read 2).
student_read_2_length = 225
check_tcr_sequencing_reach(student_read_2_length, tcr_chain='beta')
