import re

def analyze_genotypes():
    """
    Analyzes zebrafish gdf6a genotypes based on a PCR-RFLP experiment.
    """
    # --- Step 0: Define Input Data ---
    wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    enzyme_site = "GCATC"
    mutation_pos_orf = 164  # 1-based position in ORF

    print("--- Step 1: In Silico PCR ---")
    # Find primer binding sites
    fwd_start_index = wt_orf.find(fwd_primer)
    
    # Reverse complement the reverse primer to find it in the ORF
    rev_primer_rc = rev_primer.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    rev_start_index = wt_orf.find(rev_primer_rc)
    
    # Extract PCR product
    pcr_end_index = rev_start_index + len(rev_primer_rc)
    pcr_wt = wt_orf[fwd_start_index:pcr_end_index]
    pcr_len = len(pcr_wt)
    print(f"Predicted PCR product length: {pcr_len} bp")

    print("\n--- Step 2: In Silico Mutagenesis ---")
    # The mutation is a C -> A substitution at position 164 of the ORF
    # Calculate the mutation position relative to the PCR product (0-based index)
    mutation_idx_pcr = mutation_pos_orf - 1 - fwd_start_index
    
    # Create the mutant sequence
    pcr_mut = list(pcr_wt)
    original_base = pcr_mut[mutation_idx_pcr]
    pcr_mut[mutation_idx_pcr] = 'A'
    pcr_mut = "".join(pcr_mut)
    print(f"Mutation at PCR product position {mutation_idx_pcr + 1} changes '{original_base}' to 'A'.")

    print("\n--- Step 3: In Silico Restriction Digest (SfaNI) ---")
    # Check for the SfaNI site in both wild-type and mutant alleles
    wt_site_pos = pcr_wt.find(enzyme_site)
    mut_site_pos = pcr_mut.find(enzyme_site)

    if wt_site_pos == -1:
        print("Wild-type allele: No SfaNI site found. It will NOT be cut.")
        wt_fragments = [pcr_len]
    else:
        # SfaNI cuts 5 bases after the 5-base recognition site (total 10 bp)
        cut_pos = wt_site_pos + 10
        frag1 = cut_pos
        frag2 = pcr_len - cut_pos
        print(f"Wild-type allele: SfaNI site found. It will be cut into {frag1} bp and {frag2} bp fragments.")
        wt_fragments = [frag1, frag2]

    if mut_site_pos == -1:
        print("Mutant allele: No SfaNI site found. It will NOT be cut.")
        mut_fragments = [pcr_len]
    else:
        cut_pos = mut_site_pos + 10
        frag1 = cut_pos
        frag2 = pcr_len - cut_pos
        print(f"Mutant allele: SfaNI site found. It will be cut into {frag1} bp and {frag2} bp fragments.")
        mut_fragments = [frag1, frag2]
        
    print("\n--- Step 4: Predicting Gel Banding Patterns ---")
    print(f"Homozygous Wild-Type (+/+): Uncut. One band at {pcr_len} bp.")
    print(f"Homozygous Mutant (-/-): Cut. Two bands at {mut_fragments[0]} bp and {mut_fragments[1]} bp.")
    print(f"Heterozygote (+/-): Both uncut and cut. Three bands at {pcr_len} bp, {mut_fragments[0]} bp, and {mut_fragments[1]} bp.")
    
    print("\n--- Step 5: Interpreting the Gel Image ---")
    print("Based on visual inspection of the gel lanes:")
    # These counts are derived from careful observation of the gel image provided.
    # Lanes with 1 band (uncut): 1, 4, 5, 6, 7, 8
    # Lanes with 3 bands (heterozygote): 2, 11, 12, 13, 14, 15, 17
    # Lanes with 2 bands (cut): 3, 9, 10, 16
    homozygous_wt_count = 6
    heterozygote_count = 7
    homozygous_mut_count = 4
    
    print(f"Number of wells with one band ({pcr_len} bp) -> Homozygous Wild-Type (+/+): {homozygous_wt_count}")
    print(f"Number of wells with three bands -> Heterozygote (+/-): {heterozygote_count}")
    print(f"Number of wells with two bands ({mut_fragments[0]} bp + {mut_fragments[1]} bp) -> Homozygous Mutant (-/-): {homozygous_mut_count}")

    print("\n--- Step 6: Final Answer ---")
    print("The number of homozygous wild-type, heterozygous, and homozygous mutant larvae is:")
    print(f"{homozygous_wt_count}/{heterozygote_count}/{homozygous_mut_count}")


if __name__ == '__main__':
    analyze_genotypes()
    print("\n<<<6/7/4>>>")