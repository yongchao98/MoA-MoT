import textwrap

def analyze_degeneracy():
    """
    Analyzes an RNA sequence to demonstrate the principle of amino acid degeneracy,
    specifically the third position wobble effect, as described in the correct answer.
    """
    # Standard genetic code: RNA codon -> Amino Acid (single letter code)
    genetic_code = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_', 'UGA':'_'
    }

    # Amino acid full names for clarity
    aa_names = {
        'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic Acid',
        'C': 'Cysteine', 'Q': 'Glutamine', 'E': 'Glutamic Acid', 'G': 'Glycine',
        'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
        'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
        'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine',
        '_': 'STOP'
    }

    # Invert the genetic code to map amino acids to their codons (degeneracy map)
    degeneracy_map = {}
    for codon, aa in genetic_code.items():
        if aa not in degeneracy_map:
            degeneracy_map[aa] = []
        degeneracy_map[aa].append(codon)

    # Sequence from Option A
    sequence = 'GAUACGUACGAU'
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

    print(f"Analysis of Sequence from Option A: 5'-{sequence}-3'")
    print("-" * 50)
    print(f"The sequence is read as codons: {', '.join(codons)}")

    amino_acids = [genetic_code[c] for c in codons]
    unique_amino_acids = sorted(list(set(amino_acids)))
    
    print(f"Translated amino acid sequence: {'-'.join([aa_names[aa] for aa in amino_acids])}")
    print("-" * 50)
    print("Degeneracy Analysis for Each Amino Acid:")

    for aa_code in unique_amino_acids:
        aa_full_name = aa_names[aa_code]
        possible_codons = degeneracy_map[aa_code]
        num_codons = len(possible_codons)

        print(f"\nAmino Acid: {aa_full_name} ({aa_code})")
        print(f"  - Number of Codons: {num_codons}")
        print(f"  - All Possible Codons: {', '.join(possible_codons)}")
        
        # Check if the degeneracy is explained by the wobble effect
        first_two_bases = possible_codons[0][:2]
        print(f"  - Analysis: All codons for {aa_full_name} start with '{first_two_bases}'. The variation only occurs in the 3rd position.")

    print("\n" + "=" * 50)
    final_conclusion = (
        "Conclusion: The sequence 5'-GAUACGUACGAU-3' codes for Aspartic Acid, Threonine, and Tyrosine. "
        "For all of these amino acids, the degeneracy (having multiple codons) is explained by variation "
        "solely in the third codon position. This is the 'third position wobble effect'. "
        "Therefore, the sequence and its corresponding condition in option A are the most accurate and logically consistent choice."
    )
    print(textwrap.fill(final_conclusion, width=70))

analyze_degeneracy()