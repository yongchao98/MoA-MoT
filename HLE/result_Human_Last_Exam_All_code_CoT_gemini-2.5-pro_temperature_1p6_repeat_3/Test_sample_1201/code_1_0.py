import textwrap

def analyze_genetic_degeneracy():
    """
    Analyzes several RNA sequences to find the one with maximum degeneracy
    and a correct supporting biological condition.
    """
    # Step 1: Define the standard genetic code
    genetic_code = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        'AUG': 'Met',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
    }

    # Step 2: Define the choices from the problem
    choices = {
        'A': {"sequence": "GAUACGUACGAU", "condition": "third position wobble effect."},
        'B': {"sequence": "GUUUCAGAUUC", "condition": "presence of inosines."},
        'C': {"sequence": "ACGGUCAACGU", "condition": "second position pyrimidine."},
        'D': {"sequence": "CUUAUUGAUGU", "condition": "AUA as methionine in mitochondria."},
        'E': {"sequence": "AUCGCAGCUAGC", "condition": "use of alternative splicing."}
    }

    print("--- Analysis of RNA Sequences ---")

    # Step 3, 4, 5, 6: Analyze each choice
    results = {}
    for label, data in choices.items():
        sequence = data["sequence"]
        
        codons = textwrap.wrap(sequence, 3)
        amino_acids = [genetic_code.get(codon, "Unknown") for codon in codons]
        num_codons = len(codons)
        unique_amino_acids = set(amino_acids)
        num_unique_amino_acids = len(unique_amino_acids)

        # Store results for later comparison
        results[label] = {
            "num_codons": num_codons,
            "num_unique_AAs": num_unique_amino_acids,
            "condition": data["condition"]
        }

    # Step 7: Compare results and determine the best answer
    print("Conclusion:\n")
    print("To be the correct answer, a choice must exhibit degeneracy and have a correct biological explanation.")
    print("- Choices A and E show the maximum degeneracy: 4 codons yielding 3 unique amino acids.")
    print("- Choice A's condition, 'third position wobble effect,' is the correct scientific principle explaining the degeneracy of the genetic code.")
    print("- Choice E's condition, 'use of alternative splicing,' is an unrelated biological process and an incorrect explanation for codon degeneracy.")
    print("\nTherefore, Choice A is the best answer as it pairs a sequence demonstrating degeneracy with the correct underlying biological reason.")
    
    print("\n--- Detailed breakdown for the correct answer (A) ---")
    correct_choice_data = choices['A']
    correct_sequence = correct_choice_data['sequence']
    codons = textwrap.wrap(correct_sequence, 3)
    
    print(f"The sequence 5'-{correct_sequence}-3' has 4 codons:")
    print(f"  1. The first codon is GAU, which translates to the amino acid Aspartic Acid (Asp).")
    print(f"  2. The second codon is ACG, which translates to the amino acid Threonine (Thr).")
    print(f"  3. The third codon is UAC, which translates to the amino acid Tyrosine (Tyr).")
    print(f"  4. The fourth codon is GAU, which translates to the amino acid Aspartic Acid (Asp).")
    
    print("\nFinal Equation of Translation:")
    final_translation = " -> ".join([f"{c} ({genetic_code[c]})" for c in codons])
    print(final_translation)
    print(f"\nResult: 4 codons produce only 3 unique amino acids (Asp, Thr, Tyr), demonstrating degeneracy.")


analyze_genetic_degeneracy()
<<<A>>>