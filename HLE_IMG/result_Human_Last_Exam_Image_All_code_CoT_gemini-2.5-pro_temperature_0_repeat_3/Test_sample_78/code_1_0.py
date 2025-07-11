import sys
import io

# Redirect stdout to capture print statements
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_genotypes():
    """
    Performs in silico RFLP analysis to predict banding patterns and provides the final count
    based on interpreting the gel image.
    """
    # --- Provided Data ---
    # The problem states a C->A mutation at position 164. Analysis of the provided sequence
    # shows an 'A' at this position and a SfaNI site ('GCATC') ending at this position.
    # This implies the provided sequence is the MUTANT, and the wild-type (WT) has a 'C',
    # which completes the SfaNI site.
    # Therefore: WT is cut, Mutant is uncut.
    
    wt_orf_list = list("ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG")
    # Correcting the sequence to be Wild-Type (C at pos 164)
    wt_orf_list[163] = 'C'
    ORF_SEQ = "".join(wt_orf_list)
    
    FWD_PRIMER = "TTTTACGCGCTCTTCGTTTT"
    REV_PRIMER = "TTTTCCCTTGTCCACGAAAC"
    
    # --- Step 1: In Silico PCR ---
    fwd_start_index = ORF_SEQ.find(FWD_PRIMER)
    rev_comp_primer = str.maketrans('ATCG', 'TAGC')
    rev_primer_comp_seq = REV_PRIMER.translate(rev_comp_primer)[::-1]
    rev_start_index = ORF_SEQ.find(rev_primer_comp_seq)
    
    pcr_product_start = fwd_start_index
    pcr_product_end = rev_start_index + len(rev_primer_comp_seq)
    pcr_product_seq = ORF_SEQ[pcr_product_start:pcr_product_end]
    pcr_product_size = len(pcr_product_seq)

    # --- Step 2: In Silico Restriction Digest ---
    # SfaNI cuts after the C in 'GCATC'. The cut site is at position 164 of the ORF.
    cut_pos_in_orf = 164
    # Position of the cut relative to the start of the PCR product
    cut_pos_in_pcr = cut_pos_in_orf - (fwd_start_index + 1) + 1
    
    fragment1_size = cut_pos_in_pcr
    fragment2_size = pcr_product_size - fragment1_size

    print("Predicted Gel Banding Patterns:")
    print("="*35)
    
    # Homozygous Wild-Type (WT/WT)
    print("1. Homozygous Wild-Type (+/+):")
    print(f"   - The wild-type PCR product is {pcr_product_size} bp.")
    print(f"   - SfaNI cuts this product into two fragments.")
    print(f"   - Expected bands: {fragment1_size} bp and {fragment2_size} bp.")
    print("   - Gel Appearance: Two low bands, very close together (below 250 bp).")
    print("-" * 35)

    # Homozygous Mutant (mut/mut)
    print("2. Homozygous Mutant (-/-):")
    print(f"   - The C->A mutation destroys the SfaNI site.")
    print(f"   - The PCR product is not cut.")
    print(f"   - Expected band: {pcr_product_size} bp (uncut).")
    print("   - Gel Appearance: One high band (just above 250 bp).")
    print("-" * 35)

    # Heterozygote (WT/mut)
    print("3. Heterozygote (+/-):")
    print(f"   - Contains one wild-type allele and one mutant allele.")
    print(f"   - Expected bands: {pcr_product_size} bp (from mutant) and {fragment1_size} bp, {fragment2_size} bp (from wild-type).")
    print("   - Gel Appearance: Three bands (one high, two low).")
    print("="*35)

    # --- Step 3: Gel Image Analysis ---
    # By comparing the predicted patterns to the gel image, we count the lanes.
    # There are 16 sample lanes. Lane 10 appears to be a failed reaction and is excluded.
    num_homozygous_wt = 3   # Lanes with two low bands
    num_heterozygous = 10   # Lanes with three bands
    num_homozygous_mut = 2  # Lanes with one high band

    print("\nAnalysis of Gel Image (15 successful samples):")
    print(f" - Number of homozygous wild-type larvae = {num_homozygous_wt}")
    print(f" - Number of heterozygous larvae = {num_heterozygous}")
    print(f" - Number of homozygous mutant larvae = {num_homozygous_mut}")
    print("\nFinal Answer Format (WT/Het/Mut):")
    print(f"{num_homozygous_wt}/{num_heterozygous}/{num_homozygous_mut}")

# Run the analysis
analyze_genotypes()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())

# Final answer in the required format
print("<<<3/10/2>>>")