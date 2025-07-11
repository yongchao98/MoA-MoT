# The forward DNA sequence in 5' to 3' direction
dna_sequence = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"

# A dictionary to map DNA codons to single-letter amino acid codes
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', # '_' represents a Stop codon
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

# Find the starting position of the first start codon 'ATG'
start_index = dna_sequence.find('ATG')

protein_sequence = ""
# Proceed only if a start codon was found
if start_index != -1:
    # A list to store the resulting amino acids
    amino_acids = []
    
    # Iterate through the DNA sequence in steps of 3 (for each codon)
    # starting from the identified start codon
    for i in range(start_index, len(dna_sequence), 3):
        # Ensure we have a full codon of 3 bases
        if i + 3 <= len(dna_sequence):
            codon = dna_sequence[i:i+3]
            amino_acid = codon_table.get(codon, '') # Get the amino acid from the table

            # If the codon is a stop codon, terminate translation
            if amino_acid == '_':
                break
            
            # Add the translated amino acid to our list
            if amino_acid:
                 amino_acids.append(amino_acid)
        else:
            # Reached the end of the sequence before a full codon
            break
            
    # Join the list of amino acids into the final protein string
    protein_sequence = "".join(amino_acids)

# Print the final amino acid sequence
if protein_sequence:
    print(protein_sequence)
else:
    print("No start codon 'ATG' found; no protein could be translated.")
