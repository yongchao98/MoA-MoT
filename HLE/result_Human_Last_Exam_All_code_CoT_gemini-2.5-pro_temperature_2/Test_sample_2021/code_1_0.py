import math

# --- Step 1: Calculate the frequency of resistant offspring ---
self_pollination_rate = 0.05
cross_pollination_rate = 0.95
resistance_from_selfing = 0.75  # 50% wt/insert105 + 25% insert105/insert105
resistance_from_crossing = 0.50  # 50% wt/insert105

total_resistant_offspring_freq = (self_pollination_rate * resistance_from_selfing) + (cross_pollination_rate * resistance_from_crossing)
total_resistant_offspring_percent = total_resistant_offspring_freq * 100

print("### Analysis of Offspring Resistance ###")
print(f"Contribution from self-pollination: {self_pollination_rate} * {resistance_from_selfing*100}% = {self_pollination_rate * resistance_from_selfing * 100:.2f}%")
print(f"Contribution from cross-pollination: {cross_pollination_rate} * {resistance_from_crossing*100}% = {cross_pollination_rate * resistance_from_crossing * 100:.2f}%")
print(f"Total theoretical frequency of resistant offspring: {total_resistant_offspring_percent:.2f}%")
print("-" * 20)

# --- Step 2 & 3: Analyze protein activity and interaction ---
print("### Analysis of Protein Function ###")
print("E3ub-wt Activity: Co-expression with Par22 leads to a decrease in Par22 levels (from 700 to 200 units).")
print("Conclusion: E3ub-wt is an active E3 ubiquitin ligase that degrades Par22.")
print("\nE3ub-insert105 Activity: Co-expression with Par22 leads to an increase in Par22 levels (from 700 to 3000 units).")
print("Conclusion: E3ub-insert105 is not an active E3 ubiquitin ligase.")
print("\nE3ub-wt Interaction: Native MS shows a single 110 kDa peak (50 kDa Par22 + 60 kDa E3ub-wt).")
print("Conclusion: E3ub-wt interacts with Par22.")
print("\nE3ub-insert105 Interaction: Native MS shows separate peaks for Par22 and E3ub-insert105.")
print("Conclusion: Par22 cannot interact with E3ub-insert105.")
print("-" * 20)

# --- Step 4: Calculate the mass of the insertion ---
dna_insert = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc"
nucleotide_count = len(dna_insert)
amino_acid_count = nucleotide_count / 3

# Molecular weights of amino acids (in Daltons, monoisotopic)
# This dictionary maps codon to amino acid and its mass
aa_data = {
    'GAT': ('D', 115.02694), 'GAA': ('E', 129.04259), 'AAG': ('K', 128.09496),
    'TGT': ('C', 103.00919), 'TTA': ('L', 113.08406), 'CCA': ('P', 97.05276),
    'AAG': ('K', 128.09496), 'AAC': ('N', 114.04293), 'TGG': ('W', 186.07931),
    'ATC': ('I', 113.08406), 'AGA': ('R', 156.10111), 'CCA': ('P', 97.05276),
    'TTG': ('L', 113.08406), 'AAC': ('N', 114.04293), 'AGC': ('S', 87.03203),
    'TGA': ('*', 0), 'ACGAATGC': ('',0), # Placeholder, process by 3
    'GAT': ('D', 115.02694), 'GAC': ('D', 115.02694)
    # The below codon map is more robust
}

codon_map = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*',
    'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H',
    'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S',
    'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V',
    'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
    'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

aa_mass = {
    'A': 71.03711,  'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
    '*': 0
}

# Translate DNA and calculate mass
aa_sequence = ""
peptide_mass = 18.01056 # Mass of H2O is added back, as it's lost for each peptide bond but not at ends
for i in range(0, len(dna_insert), 3):
    codon = dna_insert[i:i+3].upper()
    if codon in codon_map:
        aa = codon_map[codon]
        aa_sequence += aa
        if aa in aa_mass:
            peptide_mass += aa_mass[aa]

print("### Analysis of Mass Increase ###")
print(f"The insertion has {nucleotide_count} nucleotides, which code for {int(amino_acid_count)} amino acids.")
print(f"Translated amino acid sequence: {aa_sequence}")
print(f"The calculated mass of this peptide insert is {peptide_mass:.2f} Da.")
print(f"This is approximately {peptide_mass/1000:.2f} kDa, which is rounded to 4.0 kDa.")
print("-" * 20)

print("\n### Final Conclusion ###")
print("Based on the analysis:")
print(f"1. The theoretical percentage of resistant offspring is 51.25%.")
print("2. Only E3ub-wt is an active ubiquitin ligase.")
print("3. Par22 cannot interact with the mutated E3ub-insert105 protein.")
print("4. The insertion increases the protein mass by approximately 4.0 kDa.")
print("\nThis matches all statements in option J.")