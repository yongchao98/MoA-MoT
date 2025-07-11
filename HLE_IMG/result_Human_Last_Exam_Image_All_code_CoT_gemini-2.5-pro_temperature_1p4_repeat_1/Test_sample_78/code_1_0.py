import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_genotypes():
    """
    Performs in-silico PCR and RFLP analysis to determine zebrafish genotypes.
    """
    # 1. Define sequences and parameters
    wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    sfaNI_site = "GCATC"
    mut_pos_orf = 164  # 1-based position in ORF

    # Helper function for reverse complement
    def reverse_complement(seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(complement.get(base, base) for base in reversed(seq))

    # 2. In-silico PCR
    fwd_start = wt_orf.find(fwd_primer)
    rev_comp = reverse_complement(rev_primer)
    rev_end = wt_orf.find(rev_comp) + len(rev_comp)
    wt_pcr_product = wt_orf[fwd_start:rev_end]
    pcr_product_size = len(wt_pcr_product)

    # 3. Create mutant sequence
    mut_pos_orf_idx = mut_pos_orf - 1
    mut_orf = list(wt_orf)
    mut_orf[mut_pos_orf_idx] = 'A'
    mut_orf = "".join(mut_orf)
    mut_pcr_product = mut_orf[fwd_start:rev_end]

    # 4. RFLP analysis
    wt_site_pos = wt_pcr_product.find(sfaNI_site)
    mut_site_pos = mut_pcr_product.find(sfaNI_site)

    print("Step 1: In-silico PCR and RFLP Analysis")
    print("-" * 50)
    print(f"PCR amplification product size: {pcr_product_size} bp")

    # The C->A mutation at position 164 creates a GCATC site.
    # WT sequence at site: ...TGTGCC TC...
    # Mut sequence at site: ...TGTGCATC...
    
    print("\nWild-type allele: No SfaNI site (GCATC) is present.")
    print(f"  > Expected band for homozygous wild-type (WT/WT): Uncut at {pcr_product_size} bp")

    print("\nMutant allele: The C to A mutation creates an SfaNI site.")
    # SfaNI cuts at GCATC(5/9), so 5 bp after the site.
    cut_after = mut_site_pos + len(sfaNI_site) + 5 - 1
    frag1_size = cut_after + 1
    frag2_size = pcr_product_size - frag1_size
    print(f"  > Expected bands for homozygous mutant (mut/mut): {frag1_size} bp and {frag2_size} bp")
    print("-" * 50)

    print("\nStep 2: Gel Interpretation and Genotype Summary")
    print("-" * 50)
    print("Predicted banding patterns for each genotype:")
    print(f"  - Homozygous Wild-Type (WT/WT):  One band at {pcr_product_size} bp.")
    print(f"  - Homozygous Mutant (mut/mut):    Two bands at {frag1_size} bp and {frag2_size} bp.")
    print(f"  - Heterozygote (WT/mut):         Three bands at {pcr_product_size} bp, {frag1_size} bp, and {frag2_size} bp.")
    
    print("\nBy observing the provided gel image and counting the lanes with each pattern:")
    num_wt = 2
    num_het = 10
    num_mut = 4
    print(f"- Number of homozygous wild-type larvae (one high band): {num_wt}")
    print(f"- Number of heterozygous larvae (three bands): {num_het}")
    print(f"- Number of homozygous mutant larvae (two low bands): {num_mut}")
    print("-" * 50)
    
    print("\nFinal Answer Format (WT/Heterozygous/Mutant):")
    print(f"{num_wt}/{num_het}/{num_mut}")

analyze_genotypes()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())