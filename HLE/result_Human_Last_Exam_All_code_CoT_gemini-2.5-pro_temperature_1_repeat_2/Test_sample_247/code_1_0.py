import textwrap

def simulate_tcr_sequencing():
    """
    Simulates TCR sequencing to demonstrate the challenge of capturing the
    CDR3 region with a standard 3' end protocol.
    """

    # 1. Define the structure and length of a hypothetical TCR beta transcript
    # These lengths are approximations for demonstration.
    len_v_region = 290
    len_cdr3_region = 45
    len_j_region = 50
    len_c_region = 450
    len_3_utr = 200
    len_poly_a = 50

    # The CDR3 region is what we want to sequence
    cdr3_sequence = "CASSLGQAYEQYF" # A real TCR beta CDR3 sequence example
    # Pad the sequence to its defined length
    cdr3_full = cdr3_sequence.center(len_cdr3_region, '-')

    # Build the full transcript string
    transcript = (
        "V" * len_v_region +
        cdr3_full +
        "J" * len_j_region +
        "C" * len_c_region +
        "U" * len_3_utr +
        "A" * len_poly_a
    )
    total_len = len(transcript)

    # 2. Define sequencing parameters
    read_length = 225 # As specified in the problem

    # 3. Simulate the original experiment (Standard 3' Sequencing)
    # The read starts from the 3' end (just before the PolyA tail) and goes 225bp.
    read_start_pos_from_3_end = len(transcript) - len_poly_a - read_length
    original_read = transcript[read_start_pos_from_3_end : read_start_pos_from_3_end + read_length]
    dist_cdr3_to_polya = len_j_region + len_c_region + len_3_utr
    
    # 4. Simulate the proposed solution (Targeted PCR Enrichment - Option D)
    # PCR primers are used to amplify a region including the CDR3.
    # The reverse primer binds in the Constant region. Let's say 50 bp into the C-region.
    # The sequencing read effectively starts just upstream of the V(D)J region.
    # We will simulate the read starting 50 bp upstream of the CDR3.
    targeted_read_start_pos = transcript.find(cdr3_full) - 50
    targeted_read = transcript[targeted_read_start_pos : targeted_read_start_pos + read_length]

    # 5. Print the results
    print("--- TCR Transcript Simulation ---")
    print(f"Total Transcript Length: {total_len} bp")
    print(f"CDR3 Sequence to find: {cdr3_sequence}")
    print(f"Distance from start of Poly(A) tail to CDR3: {dist_cdr3_to_polya} bp")
    print(f"Sequencing Read Length: {read_length} bp\n")

    print("--- Scenario 1: Original 3' End Sequencing ---")
    print(f"The 225 bp read starts from the 3' end.")
    print(f"Because the distance to CDR3 ({dist_cdr3_to_polya} bp) is greater than the read length ({read_length} bp), the read fails to cover the CDR3 region.")
    print(f"Is CDR3 sequence found in the read? {'Yes' if cdr3_sequence in original_read else 'No'}")
    print("Sequenced content (truncated):")
    print(textwrap.fill(original_read, width=80))
    print("\n" + "="*80 + "\n")

    print("--- Scenario 2: Targeted PCR Enrichment (Solution D) ---")
    print("A TCR-specific primer in the C-region is used for PCR amplification.")
    print("This makes the sequencing read start near the V(D)J/CDR3 region.")
    print("Is CDR3 sequence found in the read? {'Yes' if cdr3_sequence in targeted_read else 'No'}")
    print("Sequenced content (truncated):")
    print(textwrap.fill(targeted_read, width=80))

if __name__ == '__main__':
    simulate_tcr_sequencing()

<<<D>>>