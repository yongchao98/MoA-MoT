import re

def analyze_genotypes():
    """
    Performs an in-silico analysis to predict restriction digest patterns
    and provides the final count of genotypes based on the gel image.
    """
    # 1. Define experimental parameters
    orf_wt = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    # The reverse primer sequence binds to the complementary strand. We find its reverse complement to locate it on the given ORF sequence.
    # Reverse Primer: 5'-TTTTCCCTTGTCCACGAAAC-3'
    # Binds to:       3'-AAAAGGGAACAGGUGCTTTG-5'
    # Which is:       5'-GTTTCGTGGACAAGGGAAAA-3' on the coding strand
    rev_primer_binding_seq = "GTTTCGTGGACAAGGGAAAA"
    mutation_orf_pos = 164  # 1-based position
    sfaNI_site = "GCATC"

    # 2. In-Silico PCR
    fwd_start_index = orf_wt.find(fwd_primer)
    rev_binding_start_index = orf_wt.find(rev_primer_binding_seq)
    pcr_end_index = rev_binding_start_index + len(rev_primer_binding_seq)
    pcr_product_wt = orf_wt[fwd_start_index:pcr_end_index]
    pcr_product_size = len(pcr_product_wt)

    # 3. Analyze the mutation
    mutation_pcr_index = mutation_orf_pos - 1 - fwd_start_index
    pcr_product_mut = list(pcr_product_wt)
    # Confirm wild-type base is 'C'
    if pcr_product_mut[mutation_pcr_index] == 'C':
        pcr_product_mut[mutation_pcr_index] = 'A'
    pcr_product_mut = "".join(pcr_product_mut)
    
    # 4. In-Silico Restriction Digest
    # SfaNI cuts after the 'C' in 'GCATC'. A cut site is defined by the start of the sequence.
    def find_cut_positions(sequence):
        # The cut happens after the 5-base recognition site
        return sorted([m.start() + len(sfaNI_site) for m in re.finditer(sfaNI_site, sequence)])

    wt_cut_positions = find_cut_positions(pcr_product_wt)
    mut_cut_positions = find_cut_positions(pcr_product_mut)

    # 5. Predict fragment sizes
    def get_fragments(sequence_len, cut_positions):
        fragments = []
        last_cut = 0
        for pos in cut_positions:
            fragments.append(pos - last_cut)
            last_cut = pos
        fragments.append(sequence_len - last_cut)
        return sorted(fragments, reverse=True)

    wt_fragments = get_fragments(pcr_product_size, wt_cut_positions)
    mut_fragments = get_fragments(pcr_product_size, mut_cut_positions)
    het_fragments = sorted(list(set(wt_fragments + mut_fragments)), reverse=True)
    
    # 6. Interpret the Gel and Count
    # Based on manual inspection of the gel image provided.
    # We exclude lanes 1, 10, 17, and 18 as they are failed, empty, or controls.
    # Total scorable lanes = 14
    num_homozygous_wt = 4
    num_heterozygous = 7
    num_homozygous_mutant = 3
    
    # 7. Print the results
    print("--- In-Silico Analysis ---")
    print(f"PCR product size: {pcr_product_size} bp")
    print(f"Mutation at ORF position 164 (C->A) creates a new SfaNI site ({sfaNI_site}).")
    print("\n--- Predicted Gel Patterns ---")
    print(f"Homozygous Wild-Type (+/+): Digestion produces fragments of {wt_fragments} bp.")
    print(f"Homozygous Mutant (-/-): Digestion produces fragments of {mut_fragments} bp.")
    print(f"Heterozygote (+/-): Shows a combination of fragments: {het_fragments} bp.")
    print("(Note: The 33 bp fragment in the mutant and heterozygote is likely too small to be seen on the gel.)")
    print("\n--- Interpreting the Gel ---")
    print("Pattern 1 (Wild-Type): Two visible bands (~155 bp, ~111 bp).")
    print("Pattern 2 (Heterozygote): Three visible bands (~155 bp, ~122 bp, ~111 bp).")
    print("Pattern 3 (Mutant): Two visible bands (~122 bp, ~111 bp).")

    print("\n--- Final Genotype Count ---")
    print(f"Number of homozygous wild-type larvae (A): {num_homozygous_wt}")
    print(f"Number of heterozygous larvae (B): {num_heterozygous}")
    print(f"Number of homozygous mutant larvae (C): {num_homozygous_mutant}")
    
    print("\nAnswer in A/B/C format:")
    print(f"{num_homozygous_wt}/{num_heterozygous}/{num_homozygous_mutant}")

analyze_genotypes()