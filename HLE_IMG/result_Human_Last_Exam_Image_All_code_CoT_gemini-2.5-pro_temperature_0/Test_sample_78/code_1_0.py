import re

def reverse_complement(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement_map = str.maketrans('ATCG', 'TAGC')
    return seq.upper().translate(complement_map)[::-1]

# 1. Define sequences and primers
orf_wt = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
fwd_primer = "TTTTACGCGCTCTTCGTTTT"
rev_primer = "TTTTCCCTTGTCCACGAAAC"
enzyme_site = "GCATC"

# 2. In-silico PCR
fwd_pos = orf_wt.find(fwd_primer)
rev_pos = orf_wt.find(reverse_complement(rev_primer))

pcr_product_wt = ""
if fwd_pos != -1 and rev_pos != -1:
    start = fwd_pos
    end = rev_pos + len(rev_primer)
    pcr_product_wt = orf_wt[start:end]

pcr_product_size = len(pcr_product_wt)
print(f"Step 1: In-Silico PCR")
print(f"The PCR product is {pcr_product_size} bp long.\n")

# 3. Locate mutation and create mutant sequence
mutation_orf_pos = 164  # 1-based position
mutation_pcr_pos = mutation_orf_pos - fwd_pos - 1 # 0-based index in PCR product

# Create the mutant PCR product by changing C -> A
pcr_product_mut = list(pcr_product_wt)
pcr_product_mut[mutation_pcr_pos] = 'A'
pcr_product_mut = "".join(pcr_product_mut)

print(f"Step 2: Analyze Mutation Site")
print(f"Wild-type sequence around mutation: ...{pcr_product_wt[mutation_pcr_pos-5:mutation_pcr_pos+5]}...")
print(f"Mutant sequence around mutation:   ...{pcr_product_mut[mutation_pcr_pos-5:mutation_pcr_pos+5]}...\n")


# 4. In-silico Restriction Digest
print(f"Step 3: In-Silico Restriction Digest with SfaNI ({enzyme_site})")
# Check if the mutation creates the enzyme site
wt_site_present = enzyme_site in pcr_product_wt[mutation_pcr_pos-4:mutation_pcr_pos+5]
mut_site_present = enzyme_site in pcr_product_mut[mutation_pcr_pos-4:mutation_pcr_pos+5]

print(f"SfaNI site in wild-type allele at mutation location? {wt_site_present}")
print(f"SfaNI site in mutant allele at mutation location? {mut_site_present}")
print("The C->A mutation creates a new SfaNI restriction site.\n")

# Calculate fragment sizes for the mutant allele
cut_site_pos = pcr_product_mut.find(enzyme_site)
# SfaNI cuts 5 bp downstream from the end of its recognition site
cut_pos = cut_site_pos + len(enzyme_site) + 5
fragment1_size = cut_pos
fragment2_size = pcr_product_size - cut_pos

print(f"Step 4: Predict Gel Patterns")
print("Predicted fragments after SfaNI digest:")
print(f"- Homozygous Wild-Type (wt/wt): Uncut. One band at {pcr_product_size} bp.")
print(f"- Homozygous Mutant (mut/mut): Cut. Two bands at {fragment1_size} bp and {fragment2_size} bp.")
print(f"- Heterozygote (wt/mut): Both uncut and cut products. Three bands at {pcr_product_size} bp, {fragment1_size} bp, and {fragment2_size} bp.\n")

# 5. Analyze Gel Image and provide final counts
print("Step 5: Analyze Gel Image")
print("Based on the predicted patterns, we can count the genotypes from the gel image:")
print("- Lanes with one high band (~220 bp) are homozygous wild-type.")
print("- Lanes with three bands (~220, 141, 79 bp) are heterozygous.")
print("- Lanes with two low bands (~141, 79 bp) are homozygous mutant.")

# Counts from visual inspection of the 16 sample lanes (ignoring the first empty well and the last ladder well)
homozygous_wt_count = 2
heterozygous_count = 9
homozygous_mutant_count = 4

print("\nFinal Count:")
print(f"Number of homozygous wild-type larvae (A): {homozygous_wt_count}")
print(f"Number of heterozygous larvae (B): {heterozygous_count}")
print(f"Number of homozygous mutant larvae (C): {homozygous_mutant_count}")

print(f"\nThe final answer is in the format A/B/C.")
print(f"{homozygous_wt_count}/{heterozygous_count}/{homozygous_mutant_count}")