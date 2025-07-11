import sys

# Forcing UTF-8 encoding for compatibility
if sys.stdout.encoding != 'UTF-8':
    sys.stdout.reconfigure(encoding='utf-8')

def reverse_complement(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATGC", "TACG")
    return seq.upper().translate(complement_map)[::-1]

def main():
    """
    Simulates and explains why TCR CDR3 sequencing fails with 3'-end priming
    and succeeds with gene-specific priming.
    """
    # 1. Define a hypothetical TCR alpha transcript.
    # The CDR3 region is formed at the V-J junction. Lengths are for illustration.
    five_prime_utr = "GATTACAGATTACA" * 5  # 60 bp
    v_gene = "CAGTTGGCCTCAGGTGAGCC" * 4 # 80 bp
    cdr3_region = "TGTGCCAGCAGCCTAGGAGAC" # 21 bp
    j_gene = "GGAGACACAGACATTGAGGA" * 3 # 60 bp
    c_gene = "GATCTACAACAGATTGAGAG" * 5 # 100 bp
    three_prime_utr = "ATTAGCATTAGCATTAGC" * 5 # 80 bp
    poly_a_tail = "A" * 50

    tcr_mrna = five_prime_utr + v_gene + cdr3_region + j_gene + c_gene + three_prime_utr + poly_a_tail
    cdr3_start_from_5prime = len(five_prime_utr + v_gene)
    cdr3_end_from_5prime = cdr3_start_from_5prime + len(cdr3_region)
    print(f"✔️ A hypothetical TCR mRNA has been constructed. Total Length: {len(tcr_mrna)} bp.")
    print(f"    - The CDR3 region is '{cdr3_region}'.\n")


    # Bead oligo and sequencing parameters
    universal_oligo = "AAGCAGTGGTATCAACGCAGAGT"
    cell_label = "CGATCGAT"
    umi = "GATTACA"
    read1_len = 75
    read2_len = 225

    print("--- Scenario 1: Original Method (Poly-dT Priming) ---")
    # RT is primed at the poly(A) tail. cDNA is built from the 3' end of the mRNA.
    # The sequence read by Read 2 starts at the 3' end and goes "inward".
    distance_to_cdr3_from_3_end = len(three_prime_utr) + len(c_gene) + len(j_gene)
    print(f"Problem: The RT starts at the Poly(A) tail.")
    print(f"The number of bases from the 3' end to the start of the CDR3 region is: {distance_to_cdr3_from_3_end} bp.")
    print(f"The sequencing read (Read 2) has a length of {read2_len} bp.")

    if read2_len < distance_to_cdr3_from_3_end:
        print(f"Result: ❌ FAILURE. The {read2_len} bp read is too short to reach the CDR3 region, which is {distance_to_cdr3_from_3_end} bp away from the primer.")
    else:
        print(f"Result: ✅ SUCCESS. The read is long enough.")


    print("\n--- Scenario 2: Proposed Solution (C-Region Priming) ---")
    # A new bead oligo with a primer for the Constant (C) region is used.
    c_region_target_in_mrna = c_gene[20:40] # A sequence within the C-gene
    c_region_primer = reverse_complement(c_region_target_in_mrna)

    # RT now starts inside the C-gene.
    rt_start_pos_in_mrna = tcr_mrna.find(c_region_target_in_mrna)
    captured_mrna_part = tcr_mrna[:rt_start_pos_in_mrna]
    distance_to_cdr3_from_new_start = len(captured_mrna_part) - cdr3_end_from_5prime
    print("Solution: Use a gene-specific primer (GSP) that binds in the conserved C-region.")
    print(f"The RT now starts much closer to the 5' end.")
    print(f"The number of bases from the new priming site to the start of the CDR3 region is now: {distance_to_cdr3_from_new_start} bp.")
    print(f"The sequencing read (Read 2) still has a length of {read2_len} bp.")

    if read2_len > distance_to_cdr3_from_new_start:
        print(f"Result: ✅ SUCCESS. The {read2_len} bp read is easily long enough to sequence the CDR3, which is only {distance_to_cdr3_from_new_start} bp away.")
    else:
        print(f"Result: ❌ FAILURE. Even with GSP, read is not long enough (unlikely scenario).")

    print("\n--------------------------------------------------------------")
    print("Conclusion: By modifying the beads to include a capture oligo specific to the TCR constant region, the cDNA synthesis starts close enough to the CDR3 region for it to be sequenced with the existing 225 bp read length.")


if __name__ == '__main__':
    main()