import collections

def analyze_genetic_degeneracy():
    """
    Analyzes several RNA sequences and conditions to determine which best
    represents the principle of amino acid degeneracy.
    """
    # Standard genetic code mapping codons to amino acids (1-letter code).
    genetic_code = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UGC': 'C', 'UGU': 'C',
        'UGG': 'W', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
    }

    # Reverse map: Amino Acid -> List of all its codons
    aa_to_codons = collections.defaultdict(list)
    for codon, aa in genetic_code.items():
        aa_to_codons[aa].append(codon)

    choices = {
        'A': {"sequence": "GAUACGUACGAU", "condition": "third position wobble effect."},
        'B': {"sequence": "GUUUCAGAUUC", "condition": "presence of inosines."},
        'C': {"sequence": "ACGGUCAACGU", "condition": "second position pyrimidine."},
        'D': {"sequence": "CUUAUUGAUGU", "condition": "AUA as methionine in mitochondria."},
        'E': {"sequence": "AUCGCAGCUAGC", "condition": "use of alternative splicing."}
    }

    print("Analyzing each choice for amino acid degeneracy:\n")

    for key, data in choices.items():
        sequence = data["sequence"]
        condition = data["condition"]
        
        print(f"--- Option {key} ---")
        print(f"Sequence: 5'-{sequence}-3'")
        print(f"Condition: {condition}\n")

        # Step 1: Parse the sequence into codons
        codons = [sequence[i:i+3] for i in range(0, len(sequence)) if len(sequence[i:i+3]) == 3]
        
        # Step 2: Translate codons to amino acids
        amino_acids = [genetic_code.get(c, '?') for c in codons]
        peptide = "-".join(amino_acids)

        print(f"Parsed Codons: {codons}")
        print(f"Resulting Peptide: {peptide}\n")
        
        # Step 3: Analyze for demonstrated and potential degeneracy
        print("Analysis:")
        
        codon_per_aa = collections.defaultdict(list)
        for i, aa in enumerate(amino_acids):
            codon_per_aa[aa].append(codons[i])
        
        is_demonstrated = False
        for aa, c_list in codon_per_aa.items():
            if len(set(c_list)) > 1:
                print(f"  - Demonstrated Degeneracy: YES. Amino acid '{aa}' is coded by different codons: {sorted(list(set(c_list)))}.")
                is_demonstrated = True
        
        if not is_demonstrated:
            print("  - Demonstrated Degeneracy: NO. No amino acid is coded by more than one type of codon in this sequence.")

        unique_aas = sorted(list(set(amino_acids)))
        print("  - Potential Degeneracy of Produced Amino Acids:")
        for aa in unique_aas:
            if aa in aa_to_codons:
                all_codons = sorted(aa_to_codons[aa])
                print(f"    - Amino Acid '{aa}': Is {len(all_codons)}-fold degenerate. Codons: {all_codons}.")
        
        print("\n" + "="*50 + "\n")
        
    print("Final Conclusion:")
    print("The task requires finding the sequence and condition that best represent amino acid degeneracy.")
    print("- Choice E has a sequence that demonstrates degeneracy (Alanine from GCA and GCU), but the condition ('alternative splicing') is incorrect for explaining it.")
    print("- Choice A's sequence (GAU ACG UAC GAU) produces the amino acids Asp, Thr, and Tyr. While it doesn't show two different codons for one amino acid, the provided condition, 'third position wobble effect', correctly and accurately describes the biological mechanism responsible for the degeneracy of all three amino acids produced (Asp: GAU/GAC; Thr: ACU/ACC/ACA/ACG; Tyr: UAU/UAC).")
    print("\nBecause the condition in option A is a perfect match for the types of amino acids produced by its sequence, it represents the most scientifically coherent choice overall.")


if __name__ == '__main__':
    analyze_genetic_degeneracy()
<<<A>>>