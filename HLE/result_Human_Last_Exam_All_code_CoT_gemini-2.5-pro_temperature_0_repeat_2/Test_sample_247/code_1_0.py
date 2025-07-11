def simulate_tcr_sequencing():
    """
    This script simulates TCR sequencing to demonstrate why 3'-end capture fails
    to sequence the CDR3 region and how targeted priming solves the problem.
    """

    # --- 1. Define a mock TCR transcript structure and sequencing parameters ---
    # All lengths are in base pairs (bp).
    # Structure: 5' - [V-region(contains CDR3)] - [C-region] - [3' UTR] - [PolyA] - 3'
    V_REGION_LEN = 350
    C_REGION_LEN = 450
    UTR_3_LEN = 200
    
    # The CDR3 is at the 3' end of the V-region
    CDR3_LEN = 45
    CDR3_START = V_REGION_LEN - CDR3_LEN
    CDR3_END = V_REGION_LEN
    
    # Calculate total transcript length (excluding PolyA tail)
    TRANSCRIPT_LEN = V_REGION_LEN + C_REGION_LEN + UTR_3_LEN
    
    # Sequencing read length
    READ_LENGTH = 225

    print("--- TCR Transcript & Sequencing Parameters ---")
    print(f"Total Transcript Length: {TRANSCRIPT_LEN} bp")
    print(f"CDR3 Region Position: Bases {CDR3_START} to {CDR3_END}")
    print(f"Sequencing Read Length: {READ_LENGTH} bp\n")

    # --- 2. Scenario 1: Student's original method (3'-end priming) ---
    print("--- Scenario 1: Original Method (Poly-dT Priming) ---")
    # The read starts from the 3' end of the transcript and extends 225 bp "inward".
    # We calculate its coordinates relative to the 5' end of the transcript.
    read1_end_pos = TRANSCRIPT_LEN
    read1_start_pos = TRANSCRIPT_LEN - READ_LENGTH
    
    print(f"The sequencing read starts from the 3' end.")
    print(f"Equation: {read1_end_pos} (Transcript End) - {READ_LENGTH} (Read Length) = {read1_start_pos} (Read Start)")
    print(f"Resulting read covers transcript bases: {read1_start_pos} to {read1_end_pos}")

    # Check if the read covers the CDR3 region
    # Overlap exists if read_start <= cdr3_end and read_end >= cdr3_start
    if read1_start_pos <= CDR3_END and read1_end_pos >= CDR3_START:
        print("Conclusion: SUCCESS! The CDR3 region was sequenced.\n")
    else:
        print(f"Conclusion: FAILURE. The read (bases {read1_start_pos}-{read1_end_pos}) does not reach the CDR3 region (bases {CDR3_START}-{CDR3_END}).\n")

    # --- 3. Scenario 2: Proposed solution (TCR-specific priming) ---
    print("--- Scenario 2: Proposed Solution (Targeted C-Region Priming) ---")
    # A new primer is added that binds in the constant region, e.g., 50 bp into it.
    # The primer position is relative to the 5' end.
    C_REGION_START = V_REGION_LEN
    primer_pos = C_REGION_START + 50
    
    # The read starts from this primer site and extends 225 bp towards the 5' end.
    read2_end_pos = primer_pos
    read2_start_pos = primer_pos - READ_LENGTH

    print(f"A new primer binds in the Constant Region at base {primer_pos}.")
    print(f"Equation: {read2_end_pos} (Primer Position) - {READ_LENGTH} (Read Length) = {read2_start_pos} (Read Start)")
    print(f"Resulting read covers transcript bases: {read2_start_pos} to {read2_end_pos}")

    # Check if the read covers the CDR3 region
    if read2_start_pos <= CDR3_END and read2_end_pos >= CDR3_START:
        print(f"Conclusion: SUCCESS! The read (bases {read2_start_pos}-{read2_end_pos}) now covers the CDR3 region (bases {CDR3_START}-{CDR3_END}).")
    else:
        print("Conclusion: FAILURE. The CDR3 region was still not sequenced.")

# Execute the simulation
simulate_tcr_sequencing()