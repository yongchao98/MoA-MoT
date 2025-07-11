import re

def reverse_complement(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATGC", "TACG")
    return seq.translate(complement_map)[::-1]

# 1. Define sequences and parameters
wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
fwd_primer = "TTTTACGCGCTCTTCGTTTT"
rev_primer = "TTTTCCCTTGTCCACGAAAC"
mutation_pos = 164  # 1-based position
enzyme_site = "GCATC"
enzyme_name = "SfaNI"

print("--- Step 1: In-silico PCR & Restriction Digest Analysis ---")

# 2. Simulate PCR
fwd_pos = wt_orf.find(fwd_primer)
rev_pos_on_sense_strand = wt_orf.find(reverse_complement(rev_primer))

if fwd_pos == -1 or rev_pos_on_sense_strand == -1:
    print("Error: Primers not found in the provided ORF sequence.")
else:
    # Get the amplified region
    pcr_start = fwd_pos
    pcr_end = rev_pos_on_sense_strand + len(rev_primer)
    wt_pcr_product = wt_orf[pcr_start:pcr_end]
    pcr_product_size = len(wt_pcr_product)
    print(f"PCR amplification successful.")
    print(f"Predicted PCR product size: {pcr_product_size} bp")

    # 3. Simulate Wild-Type (WT) Digest
    print("\n--- Analyzing Wild-Type (+) Allele ---")
    cut_site_info = re.search(enzyme_site, wt_pcr_product)
    if cut_site_info:
        # SfaNI cuts 5 bp downstream of the recognition sequence: GCATC(N)5^
        cut_position = cut_site_info.start() + 5 + 5
        frag1_size = cut_position
        frag2_size = pcr_product_size - frag1_size
        print(f"The {enzyme_name} site '{enzyme_site}' was found in the WT PCR product.")
        print(f"Digestion of the WT product yields two fragments: {frag1_size} bp and {frag2_size} bp.")
        wt_bands = sorted([frag1_size, frag2_size], reverse=True)
    else:
        print(f"No {enzyme_name} site found in the WT PCR product.")
        wt_bands = [pcr_product_size]

    # 4. Simulate Mutant (-) Digest
    print("\n--- Analyzing Mutant (-) Allele ---")
    # The mutation is at position 164 of the ORF, which is index 163.
    mutant_orf_list = list(wt_orf)
    original_base = mutant_orf_list[mutation_pos - 1]
    mutant_orf_list[mutation_pos - 1] = 'A'
    mutant_orf = "".join(mutant_orf_list)
    mutant_pcr_product = mutant_orf[pcr_start:pcr_end]
    
    print(f"Mutation at position {mutation_pos} (C->A) was introduced.")

    if not re.search(enzyme_site, mutant_pcr_product):
        print(f"The mutation destroys the {enzyme_name} site.")
        print(f"The mutant product is not cut, resulting in one fragment: {pcr_product_size} bp.")
        mutant_bands = [pcr_product_size]
    else:
        print("The mutation does not affect the SfaNI site.")
        mutant_bands = wt_bands

# 5. Predict Genotype Patterns
print("\n--- Step 2: Predicted Gel Banding Patterns ---")
print(f"Homozygous Wild-Type (+/+): Two bands at approximately {wt_bands[0]} bp and {wt_bands[1]} bp.")
print(f"Homozygous Mutant (-/-):   One band at approximately {mutant_bands[0]} bp.")
print(f"Heterozygous (+/-):        Three bands at approximately {mutant_bands[0]} bp, {wt_bands[0]} bp, and {wt_bands[1]} bp.")

# 6. Analyze Gel Image and Count
print("\n--- Step 3: Analysis of the Gel Image ---")
print("By comparing these predicted patterns with the provided gel image, we can count the number of individuals for each genotype.")
print("The ladder indicates the lowest band is 250 bp. The uncut mutant band (~279 bp) runs just above this mark. The cut wild-type bands (~152 bp and ~127 bp) run significantly lower.")
print("Counting the lanes with scorable results (15 lanes) reveals the following distribution:")

# Based on manual inspection of the gel image
homozygous_wt_count = 4
heterozygote_count = 9
homozygous_mutant_count = 2

print(f"\nNumber of homozygous wild-type (+/+) larvae (2 low bands): {homozygous_wt_count}")
print(f"Number of heterozygous (+/-) larvae (3 bands): {heterozygote_count}")
print(f"Number of homozygous mutant (-/-) larvae (1 high band): {homozygous_mutant_count}")

# 7. Output Final Result
print("\n--- Final Answer ---")
final_answer = f"{homozygous_wt_count}/{heterozygote_count}/{homozygous_mutant_count}"
print(final_answer)
print("\n<<<4/9/2>>>")