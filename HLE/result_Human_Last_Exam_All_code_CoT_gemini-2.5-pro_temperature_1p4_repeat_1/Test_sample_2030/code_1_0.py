def translate_dna(seq):
    """Translates a DNA sequence into a peptide sequence."""
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GUG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    peptide = ""
    for i in range(0, len(seq) - (len(seq) % 3), 3):
        codon = seq[i:i+3]
        peptide += codon_table.get(codon, 'X')
    return peptide

def get_reverse_complement(dna_seq):
    """Returns the reverse complement of a DNA sequence."""
    complement_map = str.maketrans('ATCG', 'TAGC')
    complement = dna_seq.translate(complement_map)
    return complement[::-1]

# Step 1: Define original sequence and get reverse complement
original_seq_raw = "CTT CCC CGC ACA AGT GGT"
original_seq = original_seq_raw.replace(" ", "")
rev_comp_seq = get_reverse_complement(original_seq)

print("--- Step 1: Translating all 6 Reading Frames ---")
print(f"Original Sequence: 5' {original_seq_raw} 3'")
print(f"Reverse Complement: 5' {rev_comp_seq[:3]} {rev_comp_seq[3:6]} {rev_comp_seq[6:9]} {rev_comp_seq[9:12]} {rev_comp_seq[12:15]} {rev_comp_seq[15:18]} 3'")
print("-" * 20)

frames = {}
frames['F1'] = (original_seq, translate_dna(original_seq))
frames['F2'] = (original_seq[1:], translate_dna(original_seq[1:]))
frames['F3'] = (original_seq[2:], translate_dna(original_seq[2:]))
frames['F4'] = (rev_comp_seq, translate_dna(rev_comp_seq))
frames['F5'] = (rev_comp_seq[1:], translate_dna(rev_comp_seq[1:]))
frames['F6'] = (rev_comp_seq[2:], translate_dna(rev_comp_seq[2:]))

for i in range(1, 7):
    name = f'F{i}'
    print(f"Frame {i}: {frames[name][1]}")

print("\n--- Step 2: Identifying the Target Frame ---")
all_peptides = "".join([frames[name][1] for name in frames])
unique_amino_acids_per_frame = {}

for name, (dna, peptide) in frames.items():
    other_peptides = "".join([p for n, (d, p) in frames.items() if n != name])
    uniques = [aa for aa in set(peptide) if aa not in other_peptides and aa != '*']
    unique_amino_acids_per_frame[name] = uniques

target_frame_name = None
for name, uniques in unique_amino_acids_per_frame.items():
    if len(uniques) == 2:
        target_frame_name = name
        break

target_dna_seq = frames[target_frame_name][0]
target_peptide = frames[target_frame_name][1]
unique_AAs = unique_amino_acids_per_frame[target_frame_name]

print(f"Found Frame {target_frame_name} with 2 unique amino acids: {', '.join(unique_AAs)}")
print(f"Peptide sequence: {target_peptide}")

print("\n--- Step 3: Applying SNP Rules ---")
# The unique amino acids are Phenylalanine (F) and Tryptophan (W).
# F (non-polar) must be changed to Cysteine (Cys).
# W (considered polar for this problem) must be changed to a Stop codon.

original_codons = [target_dna_seq[i:i+3] for i in range(0, len(target_dna_seq) - (len(target_dna_seq)%3), 3)]
print(f"Original Codons: {' '.join(original_codons)}")

# Finding the changes
modified_codons = list(original_codons)
# Change F (TTC) -> Cys (TGC)
modified_codons[0] = "TGC"
# Change W (TGG) -> Stop (TGA)
modified_codons[4] = "TGA"

print(f"Original F codon: TTC -> Modified Cys codon: {modified_codons[0]}")
print(f"Original W codon: TGG -> Modified Stop codon: {modified_codons[4]}")
print(f"Modified peptide: {translate_dna(''.join(modified_codons))}")


print("\n--- Step 4 & 5: Constructing Modified Sequence and Designing Oligo ---")
# The oligo should bind to the sequence translated into amino acids, which is everything before the stop codon.
oligo_target_sequence = "".join(modified_codons[:-1])
print(f"Modified DNA sequence for oligo design (pre-stop codon): 5' {oligo_target_sequence} 3'")

# The oligo sequence is the reverse complement of the target.
oligo_sequence = get_reverse_complement(oligo_target_sequence)

print("\n--- Final Answer ---")
print("The final oligo binds to the modified sequence before the stop codon.")
print(f"The sequence for this oligo is: 5' {oligo_sequence} 3'")