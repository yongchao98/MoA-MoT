import collections

def analyze_degeneracy():
    """
    Analyzes RNA sequences for amino acid degeneracy and evaluates the correctness
    of paired biological conditions.
    """
    # Standard RNA codon table: Codon -> 1-Letter Amino Acid Code
    RNA_CODON_TABLE = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', # Isoleucine
        'AUG':'M', # Methionine
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T', # Threonine
        'AAC':'N', 'AAU':'N', # Asparagine
        'AAA':'K', 'AAG':'K', # Lysine
        'AGC':'S', 'AGU':'S', # Serine
        'AGA':'R', 'AGG':'R', # Arginine
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L', # Leucine
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P', # Proline
        'CAC':'H', 'CAU':'H', # Histidine
        'CAA':'Q', 'CAG':'Q', # Glutamine
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R', # Arginine
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V', # Valine
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A', # Alanine
        'GAC':'D', 'GAU':'D', # Aspartic Acid
        'GAA':'E', 'GAG':'E', # Glutamic Acid
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G', # Glycine
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S', # Serine
        'UUC':'F', 'UUU':'F', # Phenylalanine
        'UUA':'L', 'UUG':'L', # Leucine
        'UAC':'Y', 'UAU':'Y', # Tyrosine
        'UGC':'C', 'UGU':'C', # Cysteine
        'UGG':'W', # Tryptophan
        'UAA':'*', 'UAG':'*', 'UGA':'*' # Stop
    }

    # Calculate how many codons code for each amino acid (degeneracy level)
    AA_DEGENERACY = collections.Counter(RNA_CODON_TABLE.values())

    # Data from the answer choices
    options = {
        'A': {'sequence': 'GAUACGUACGAU', 'condition': 'third position wobble effect.'},
        'B': {'sequence': 'GUUUCAGAUUC', 'condition': 'presence of inosines.'},
        'C': {'sequence': 'ACGGUCAACGU', 'condition': 'second position pyrimidine.'},
        'D': {'sequence': 'CUUAUUGAUGU', 'condition': 'AUA as methionine in mitochondria.'},
        'E': {'sequence': 'AUCGCAGCUAGC', 'condition': 'use of alternative splicing.'}
    }

    print("Analyzing each option for sequence degeneracy and condition validity...\n")

    for option_letter, data in options.items():
        sequence = data['sequence']
        condition = data['condition']
        
        # Split sequence into codons, ignoring trailing nucleotides
        codons = [sequence[i:i+3] for i in range(0, len(sequence) - len(sequence) % 3, 3)]
        
        # Translate codons to amino acids
        amino_acids = [RNA_CODON_TABLE.get(c, '?') for c in codons]
        
        # Perform analysis
        num_codons = len(codons)
        unique_amino_acids = set(amino_acids)
        num_unique_aa = len(unique_amino_acids)
        
        # Calculate average degeneracy score for the translated amino acids
        degeneracy_scores = [AA_DEGENERACY.get(aa, 0) for aa in amino_acids]
        avg_degeneracy_score = sum(degeneracy_scores) / len(degeneracy_scores) if degeneracy_scores else 0

        # Print the detailed analysis for the current option
        print(f"--- Option {option_letter} ---")
        print(f"Sequence: 5'-{sequence}-3'")
        print(f"Codons: {' '.join(codons)}")
        print(f"Resulting Amino Acids: {'-'.join(amino_acids)}")
        print(f"Total Codons = {num_codons}")
        print(f"Unique Amino Acids = {num_unique_aa}")
        print(f"Average Degeneracy of Constituent Amino Acids = {avg_degeneracy_score:.2f}")
        print(f"Paired Condition: {condition}")
        
        # Evaluate the validity of the paired condition
        validity_explanation = ""
        if option_letter == 'A':
            validity_explanation = ("Condition is CORRECT. The 'wobble' in the third position of a codon is the primary mechanism for degeneracy in the genetic code.")
        elif option_letter == 'B':
            validity_explanation = ("Condition is INCORRECT. Inosine is a modified base found in tRNA anticodons, not typically in mRNA, and is not present in this sequence.")
        elif option_letter == 'C':
            validity_explanation = ("Condition is INCORRECT. This is not a recognized general principle for degeneracy. While the second base is critical, this statement is too vague and not a mechanism like wobble.")
        elif option_letter == 'D':
            validity_explanation = ("Condition is INCORRECT for this sequence. The condition refers to codon AUA in mitochondria, but the sequence contains AUU. They are different codons.")
        elif option_letter == 'E':
            validity_explanation = ("Condition is INCORRECT. Alternative splicing is a pre-mRNA processing event that creates different mRNAs from one gene; it is unrelated to the process of translating a given mRNA sequence into protein.")
        
        print(f"Evaluation: {validity_explanation}\n")
    
    # Final conclusion based on the complete analysis
    print("--- Final Conclusion ---")
    print("To be the correct answer, both the sequence and its paired condition must be valid.")
    print("Sequence E shows a high degeneracy score (4.25) and a clear example of wobble (GCA and GCU code for Alanine). However, its paired condition ('alternative splicing') is incorrect.")
    print("Sequence A shows degeneracy by producing 3 unique amino acids from 4 codons (Asp-Thr-Tyr-Asp). Crucially, its paired condition ('third position wobble effect') is the correct and fundamental scientific principle that explains this phenomenon.")
    print("\nBecause Option A is the only choice that provides a valid pair of a degenerate sequence and a correct explanatory condition, it is the correct answer.")

# Run the analysis
analyze_degeneracy()