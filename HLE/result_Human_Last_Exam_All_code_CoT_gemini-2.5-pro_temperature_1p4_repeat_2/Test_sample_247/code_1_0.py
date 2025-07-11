import sys

def simulate_tcr_sequencing():
    """
    This script simulates TCR sequencing to demonstrate why poly(dT) priming
    fails to capture the CDR3 region with a 225 bp read, while priming in
    the Constant (C) region succeeds.
    """

    # Define representative lengths (in bases) for a TCR beta chain mRNA
    LEN_5_UTR = 50
    LEN_LEADER = 60
    LEN_V_REGION = 300  # Variable region
    LEN_D_REGION = 15   # Diversity region
    LEN_J_REGION = 45   # Joining region
    LEN_C_REGION = 400  # Constant region
    LEN_3_UTR = 200
    LEN_POLYA = 100

    # The CDR3 region is formed at the V-D-J junction.
    # We calculate its position relative to the start of the PolyA tail.
    # A position of 0 is the first base of the polyA tail.
    
    # Distance from PolyA start to the end of the C-region
    end_of_c_region_pos = LEN_3_UTR
    # Distance from PolyA start to the start of the CDR3 (at the VDJ junction)
    cdr3_start_pos_from_polya = LEN_3_UTR + LEN_C_REGION + LEN_J_REGION + LEN_D_REGION
    # Distance from PolyA start to the end of the CDR3
    cdr3_end_pos_from_polya = LEN_3_UTR + LEN_C_REGION
    
    print("--- TCR mRNA Model ---")
    print(f"CDR3 region is located between base {cdr3_end_pos_from_polya} and {cdr3_start_pos_from_polya} upstream from the Poly(A) tail.")
    
    # Sequencing parameters from the problem description
    READ_2_LENGTH = 225
    print(f"Sequencing Read 2 Length: {READ_2_LENGTH} bp\n")

    # --- SCENARIO 1: Original Method (Poly-dT Priming) ---
    print("--- Scenario 1: Original Method (Poly-dT Priming) ---")
    # Reverse transcription and sequencing start from the 3' end (PolyA tail).
    # The sequenced region covers bases 0 to 225 upstream from the primer.
    farthest_base_reached_poly_dt = READ_2_LENGTH
    print(f"Sequencing from the Poly(A) tail covers bases 0 to {farthest_base_reached_poly_dt}.")
    
    # Check if the sequenced region reaches the CDR3
    if farthest_base_reached_poly_dt < cdr3_end_pos_from_polya:
        print(f"Result: FAILED. The read stops at base {farthest_base_reached_poly_dt}, which is {cdr3_end_pos_from_polya - farthest_base_reached_poly_dt} bases short of reaching the CDR3 region.\n")
    else:
        print("Result: SUCCESS. The CDR3 region would be captured.\n")

    # --- SCENARIO 2: Proposed Solution (C-Region Priming) ---
    print("--- Scenario 2: Proposed Solution (C-Region Priming as in Answer C) ---")
    # A new primer is added that binds within the Constant (C) region.
    # Let's assume it binds 50 bases into the C-region (from its 3' end).
    primer_pos_in_c = 50
    primer_pos_from_polya = LEN_3_UTR + primer_pos_in_c
    
    print(f"A new TCR-specific primer binds at base {primer_pos_from_polya} (within the Constant Region).")
    
    # The sequenced region now covers bases from the primer up to READ_2_LENGTH farther.
    farthest_base_reached_c_primer = primer_pos_from_polya + READ_2_LENGTH
    print(f"Sequencing from this primer covers the region from base {primer_pos_from_polya} to {farthest_base_reached_c_primer}.")

    # Check if the sequenced region reaches and covers the CDR3
    if farthest_base_reached_c_primer > cdr3_start_pos_from_polya:
        coverage = farthest_base_reached_c_primer - cdr3_end_pos_from_polya
        print(f"Result: SUCCESS. The read reaches past the entire CDR3 region with {coverage} bp to spare.")
        print("This ensures the full V(D)J/CDR3 region is sequenced.")
    else:
        print("Result: FAILED. Even with C-region priming, the read is not long enough in this model.")

# Run the simulation
simulate_tcr_sequencing()
