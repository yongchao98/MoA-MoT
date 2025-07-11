import sys

def solve():
    """
    This script simulates the sequencing of a T-cell receptor (TCR) transcript
    to determine if the CDR3 region can be sequenced under two different lab protocols.
    """

    # --- Define parameters for a typical T-cell receptor alpha (TRA) chain mRNA ---
    # Length of the 3' Untranslated Region (UTR) in nucleotides
    len_3_prime_utr = 200
    # Length of the Constant (C) region in nucleotides
    len_constant_region = 411
    # Approximate distance from the start of the Constant region to the CDR3 loop
    distance_from_c_region_to_cdr3 = 150
    # Illumina Read 2 length specified by the student
    read_2_length = 225

    # --- Scenario 1: Original setup with Poly(dT) capture ---
    print("--- Analysis of Original Poly(dT) Capture Method ---")
    print("This method uses a poly(dT) oligo on the bead to capture the mRNA's poly(A) tail.")
    
    # To sequence the CDR3, Read 2 must be long enough to cover the entire 3' UTR,
    # the entire Constant region, and the region between the C-region and the CDR3.
    distance_to_cdr3_from_polyA = len_3_prime_utr + len_constant_region + distance_from_c_region_to_cdr3

    print(f"The total distance from the poly(A) tail priming site to the CDR3 region is:")
    print(f"{len_3_prime_utr} bp (3' UTR) + {len_constant_region} bp (C-region) + {distance_from_c_region_to_cdr3} bp (J-region etc.) = {distance_to_cdr3_from_polyA} bp")
    print(f"The available sequencing Read 2 length is: {read_2_length} bp")

    if read_2_length >= distance_to_cdr3_from_polyA:
        print("\nResult: SUCCESS. The 225 bp read is long enough to sequence the CDR3 region.")
    else:
        print("\nResult: FAILURE. The 225 bp read is too short to reach the CDR3 region.")
        print(f"Equation: {read_2_length} bp (Read Length) < {distance_to_cdr3_from_polyA} bp (Required Length)")

    print("\n" + "="*60 + "\n")

    # --- Scenario 2: Proposed solution with TCR-specific primer (Option C) ---
    print("--- Analysis of Proposed TCR-Specific Primer Method (Option C) ---")
    print("This method modifies the bead oligo to be a specific primer that binds within the TCR Constant region.")
    
    # Assume the new primer binds 50 bp from the 3' end of the Constant region.
    # This moves the start site of reverse transcription much closer to the CDR3.
    gsp_position_from_c_region_end = 50
    distance_to_cdr3_from_gsp = (len_constant_region - gsp_position_from_c_region_end) + distance_from_c_region_to_cdr3

    print(f"The new distance from the specific priming site to the CDR3 region is:")
    print(f"({len_constant_region} - {gsp_position_from_c_region_end}) bp (Remaining C-region) + {distance_from_c_region_to_cdr3} bp (J-region etc.) = {distance_to_cdr3_from_gsp} bp")
    print(f"The available sequencing Read 2 length is: {read_2_length} bp")
    
    if read_2_length >= distance_to_cdr3_from_gsp:
        print("\nResult: SUCCESS. The 225 bp read is now long enough to sequence the CDR3 region.")
        print(f"Equation: {read_2_length} bp (Read Length) >= {distance_to_cdr3_from_gsp} bp (Required Length)")
    else:
        print("\nResult: FAILURE. Even with the new primer, the read is too short.")

solve()