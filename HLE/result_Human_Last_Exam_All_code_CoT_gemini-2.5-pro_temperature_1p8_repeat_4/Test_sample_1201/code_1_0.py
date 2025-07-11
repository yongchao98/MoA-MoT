def analyze_degeneracy():
    """
    Analyzes a given RNA sequence to demonstrate amino acid degeneracy
    based on the third position wobble effect.
    """
    # Standard genetic code: codon -> amino acid (3-letter code)
    genetic_code = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    # Reverse map: amino acid -> list of codons
    aa_to_codons = {}
    for codon, aa in genetic_code.items():
        if aa not in aa_to_codons:
            aa_to_codons[aa] = []
        aa_to_codons[aa].append(codon)

    # Sequence from choice A
    sequence = "GAUACGUACGAU"
    print(f"Analyzing sequence: 5'-{sequence}-3'\n")

    # Process the sequence into codons
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

    # Analyze each codon
    for codon in codons:
        amino_acid = genetic_code.get(codon, 'Unknown')
        if amino_acid != 'Unknown':
            all_codons_for_aa = aa_to_codons[amino_acid]
            degeneracy = len(all_codons_for_aa)
            
            print(f"Codon: {codon}")
            print(f" -> Translates to: {amino_acid}")
            print(f" -> Amino Acid Degeneracy: {degeneracy}")
            print(f" -> All codons for {amino_acid}: {', '.join(all_codons_for_aa)}")
            if degeneracy > 1:
                print("    (This demonstrates the wobble effect at the 3rd position.)")
            print("-" * 30)

    print("\nConclusion:")
    print("Choice A is the correct answer because it is the only option that pairs a sequence with the correct and fundamental biological principle explaining degeneracy: the 'third position wobble effect'. The other options provide incorrect or irrelevant conditions.")

analyze_degeneracy()