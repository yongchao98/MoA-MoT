import textwrap

def reverse_complement(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement_map.get(base, 'N') for base in reversed(seq))

def main():
    # 1. Define the components of a simplified TCR-beta mRNA transcript.
    # The CDR3 region is far from the 3' Poly-A tail.
    leader_sequence = "ATGTTTGTACAG"
    v_region = "GTGTCCAGAGGAGACAGGGCTTACTGTTCTACCAGGAGAAC"
    cdr3_region = "TGTGCCAGCAGCTTAGCGGGACAGGGGGCTACGAGCAGTACTTC"  # The important region to sequence
    j_region = "GGCACCCGGCTGACCGTGCTAG"
    c_region = "AGGACCTGAAAAACGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCAC"
    poly_a_tail = "A" * 30

    tcr_mrna = leader_sequence + v_region + cdr3_region + j_region + c_region + poly_a_tail

    print("--- Step 1: TCR mRNA structure ---")
    print(f"Total mRNA Length: {len(tcr_mrna)} bp")
    print(f"Location of CDR3 start: {tcr_mrna.find(cdr3_region)} bp from 5' end")
    print(f"Distance of CDR3 from Poly(A) tail start: {len(j_region) + len(c_region)} bp\n")
    print("This distance is too long for a standard 225 bp sequencing read from the 3' end.\n")

    # 2. Define the oligo on the bead for initial capture and the resulting cDNA.
    universal_adapter = "GATTACATATACGAGCATAC"
    umi = "ACGTACGT"
    cell_label = "AAGCTAGCTT"
    poly_t_primer = "T" * 30

    # The reverse transcriptase creates the first strand of cDNA from the mRNA template.
    full_cdna_from_mrna = reverse_complement(tcr_mrna)
    full_cdna_on_bead = universal_adapter + cell_label + umi + poly_t_primer + full_cdna_from_mrna

    print("--- Step 2: Initial cDNA Synthesis (Same for all mRNAs) ---")
    print(f"Full cDNA on bead (first 100bp): {full_cdna_on_bead[:100]}...\n")

    # 3. Design primers for targeted PCR as suggested in option D.
    # One primer binds to the universal adapter on the bead.
    # The other (reverse) primer binds to the conserved TCR constant region (C-region).
    pcr_primer_forward = universal_adapter  # Binds to the bead adapter
    pcr_primer_reverse_target = reverse_complement(c_region[:25]) # Primer binds to start of C-region on cDNA

    print("--- Step 3: Targeted PCR Amplification (The Solution) ---")
    print(f"Forward Primer (targets bead adapter): {pcr_primer_forward}")
    print(f"Reverse Primer (targets C-Region on cDNA): {pcr_primer_reverse_target}\n")

    # 4. Simulate the PCR to create the final product for sequencing.
    forward_primer_start_index = full_cdna_on_bead.find(pcr_primer_forward)
    reverse_primer_start_index = full_cdna_on_bead.find(pcr_primer_reverse_target)

    if forward_primer_start_index != -1 and reverse_primer_start_index != -1:
        # The amplicon is the sequence between the two primers.
        amplicon_start = forward_primer_start_index
        amplicon_end = reverse_primer_start_index + len(pcr_primer_reverse_target)
        tcr_amplicon = full_cdna_on_bead[amplicon_start:amplicon_end]

        # The sequenced part of the amplicon will contain the CDR3.
        # Let's check if the reverse complement of the CDR3 is in our amplicon.
        cdr3_in_cdna = reverse_complement(cdr3_region)
        is_cdr3_present = cdr3_in_cdna in tcr_amplicon

        print("--- Step 4: Final TCR Amplicon for Sequencing ---")
        print(f"Final Amplicon Length: {len(tcr_amplicon)} bp (suitable for sequencing)")
        print("Final Amplicon Sequence (shortened):")
        print(textwrap.fill(tcr_amplicon, 80))
        print(f"\nIs CDR3 region present in the final amplicon? {is_cdr3_present}")
        print("\nConclusion: By using a TCR-specific primer in the C-region for PCR, we create a library")
        print("that is enriched for TCRs and specifically includes the crucial V(D)J/CDR3 region.")

if __name__ == '__main__':
    main()
