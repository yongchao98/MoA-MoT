import collections

def analyze_degeneracy():
    """
    Analyzes the amino acid degeneracy for given RNA sequences.
    This function defines the standard genetic code, calculates the degeneracy
    for each amino acid, and then processes a list of RNA sequences to
    determine the degeneracy of the amino acids they code for.
    """
    
    # Standard genetic code: RNA codon -> Amino Acid (single-letter code)
    # '*' denotes a Stop codon.
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    # Calculate degeneracy for each amino acid
    amino_acid_counts = collections.Counter(codon_table.values())
    degeneracy_map = {codon: amino_acid_counts[aa] for codon, aa in codon_table.items()}

    # Full name mapping for amino acids
    aa_names = {
        'F': 'Phe', 'L': 'Leu', 'I': 'Ile', 'M': 'Met', 'V': 'Val',
        'S': 'Ser', 'P': 'Pro', 'T': 'Thr', 'A': 'Ala', 'Y': 'Tyr',
        'H': 'His', 'Q': 'Gln', 'N': 'Asn', 'K': 'Lys', 'D': 'Asp',
        'E': 'Glu', 'C': 'Cys', 'W': 'Trp', 'R': 'Arg', 'G': 'Gly',
        '*': 'Stop'
    }

    sequences = {
        "A": "5'-GAUACGUACGAU-3'",
        "B": "5'-GUUUCAGAUUC-3'",
        "C": "5'-ACGGUCAACGU-3'",
        "D": "5'-CUUAUUGAUGU-3'",
        "E": "5'-AUCGCAGCUAGC-3'",
    }
    
    print("Analysis of Amino Acid Degeneracy in RNA Sequences\n")

    for choice, rna_full in sequences.items():
        # Extract the raw sequence
        rna_seq = rna_full.split("'")[1]
        
        # Split into codons (of 3 bases each)
        codons = [rna_seq[i:i+3] for i in range(0, len(rna_seq), 3) if i+3 <= len(rna_seq)]
        
        print(f"Choice {choice}: {rna_full}")
        print("-" * (len(rna_full) + 9))
        
        if not codons:
            print("No full codons in sequence.\n")
            continue
            
        print(f"{'Codon':<6} | {'Amino Acid':<12} | {'Degeneracy (Codons for AA)'}")
        print("-" * 45)
        
        total_degeneracy = 0
        for codon in codons:
            amino_acid_code = codon_table.get(codon, '?')
            amino_acid_name = aa_names.get(amino_acid_code, 'Unknown')
            degeneracy = degeneracy_map.get(codon, 0)
            total_degeneracy += degeneracy
            
            # The final equation should output each number
            # Here we are outputting the degeneracy for each codon translated.
            print(f"{codon:<6} | {amino_acid_name:<12} | {degeneracy}")

        avg_degeneracy = total_degeneracy / len(codons)
        print(f"Total Degeneracy Sum = {total_degeneracy}")
        print(f"Average Degeneracy = {avg_degeneracy:.2f}\n")

analyze_degeneracy()