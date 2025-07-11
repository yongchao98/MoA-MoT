def analyze_degeneracy():
    """
    Analyzes the RNA sequence from option A to demonstrate amino acid degeneracy.
    """
    # Standard genetic code mapping amino acids to their codons
    amino_acid_codons = {
        'I': ['AUU', 'AUC', 'AUA'],
        'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
        'V': ['GUU', 'GUC', 'GUA', 'GUG'],
        'F': ['UUU', 'UUC'],
        'M': ['AUG'],
        'C': ['UGU', 'UGC'],
        'A': ['GCU', 'GCC', 'GCA', 'GCG'],
        'G': ['GGU', 'GGC', 'GGA', 'GGG'],
        'P': ['CCU', 'CCC', 'CCA', 'CCG'],
        'T': ['ACU', 'ACC', 'ACA', 'ACG'],
        'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
        'Y': ['UAU', 'UAC'],
        'W': ['UGG'],
        'Q': ['CAA', 'CAG'],
        'N': ['AAU', 'AAC'],
        'H': ['CAU', 'CAC'],
        'E': ['GAA', 'GAG'],
        'D': ['GAU', 'GAC'],
        'K': ['AAA', 'AAG'],
        'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'Stop': ['UAA', 'UAG', 'UGA']
    }

    # Reverse map for quick codon-to-amino-acid lookup
    codon_to_aa = {codon: aa for aa, codons in amino_acid_codons.items() for codon in codons}

    # Sequence from Option A
    rna_sequence = "GAUACGUACGAU"
    
    print(f"Analyzing RNA sequence: 5'-{rna_sequence}-3'")
    
    # Split sequence into codons
    codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]
    print(f"Codons: {codons}")

    translated_peptide = []
    degeneracy_values = []
    
    for codon in codons:
        amino_acid = codon_to_aa[codon]
        degeneracy = len(amino_acid_codons[amino_acid])
        
        translated_peptide.append(amino_acid)
        degeneracy_values.append(degeneracy)
        
        print(f"Codon {codon} -> Amino Acid {amino_acid} (Degeneracy: {degeneracy})")

    # The "equation" is the sum of the individual degeneracy values
    total_degeneracy_score = sum(degeneracy_values)
    
    print("\nSummary:")
    print(f"The sequence translates to the peptide: {'-'.join(translated_peptide)}")
    print("The degeneracy level for each amino acid in the sequence contributes to the overall degeneracy.")
    equation_str = " + ".join(map(str, degeneracy_values))
    print(f"Total Degeneracy Score = {equation_str} = {total_degeneracy_score}")
    print("\nConclusion: Option A is correct because it pairs a representative sequence with the correct underlying principle, the 'third position wobble effect'.")

# Execute the analysis
analyze_degeneracy()