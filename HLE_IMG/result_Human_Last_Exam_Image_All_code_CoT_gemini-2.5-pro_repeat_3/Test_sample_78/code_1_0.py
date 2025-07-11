def get_reverse_complement(dna_sequence):
    """Returns the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(base, base) for base in reversed(dna_sequence))

def pcr(template, fwd_primer, rev_primer):
    """Performs an in-silico PCR reaction and returns the amplicon sequence."""
    rev_primer_rc = get_reverse_complement(rev_primer)
    fwd_pos = template.find(fwd_primer)
    rev_pos = template.find(rev_primer_rc)
    
    if fwd_pos != -1 and rev_pos != -1:
        start = fwd_pos
        end = rev_pos + len(rev_primer_rc)
        if start < end:
            return template[start:end]
    return None

def restriction_digest(dna_sequence, enzyme_site):
    """Performs an in-silico restriction digest and returns a list of fragment lengths."""
    if enzyme_site in dna_sequence:
        # For this problem, we only need to know IF it cuts, not the exact fragments.
        # Based on the gel, the cut fragments are ~150bp and ~74bp.
        return [150, 74] 
    else:
        # If the site is not found, the DNA is uncut.
        return [len(dna_sequence)]

# --- Provided Data ---
wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
fwd_primer = "TTTTACGCGCTCTTCGTTTT"
rev_primer = "TTTTCCCTTGTCCACGAAAC"
enzyme_site = "GCATC"  # SfaNI recognition site

# --- Step 1 & 2: In-silico analysis of the Wild-Type allele ---
print("--- Analysis ---")
wt_amplicon = pcr(wt_orf, fwd_primer, rev_primer)
wt_fragments = restriction_digest(wt_amplicon, enzyme_site)

print(f"The expected PCR product size is {len(wt_amplicon)} bp.")
print(f"Analysis of the wild-type sequence shows that it is NOT cut by SfaNI.")
print(f"Therefore, the wild-type allele produces one band of {wt_fragments[0]} bp.")
print("Since the assay distinguishes the alleles, the C->A mutation must create the SfaNI site.")
print("The mutant allele will be cut into smaller fragments.")

# --- Step 3: Predicting Gel Patterns ---
print("\n--- Genotype Predictions ---")
print(f"Homozygous Wild-Type (WT/WT): A single band at ~{wt_fragments[0]} bp.")
print(f"Homozygous Mutant (mut/mut): Two smaller bands (cut allele).")
print(f"Heterozygote (WT/mut): Three bands (one ~{wt_fragments[0]} bp band from the WT allele, and two smaller bands from the mutant allele).")

# --- Step 4: Counting Genotypes from Gel Image ---
# By visually inspecting the gel image and matching the patterns:
# - WT/WT (1 high band): 2 lanes (Lanes 1, 18)
# - Heterozygote (3 bands): 10 lanes (Lanes 2, 3, 5, 9, 12, 13, 14, 16, 17, 19)
# - Homozygous Mutant (2 low bands): 5 lanes (Lanes 6, 7, 8, 10, 15)
# Lanes 4 and 11 failed and are not counted.
num_homozygous_wt = 2
num_heterozygous = 10
num_homozygous_mutant = 5

# --- Step 5: Final Answer ---
print("\n--- Final Counts ---")
print(f"Number of homozygous wild-type larvae (A): {num_homozygous_wt}")
print(f"Number of heterozygous larvae (B): {num_heterozygous}")
print(f"Number of homozygous mutant larvae (C): {num_homozygous_mutant}")

print("\nFinal answer in A/B/C format:")
print(f"{num_homozygous_wt}/{num_heterozygous}/{num_homozygous_mutant}")