def suggest_mutagenesis():
    """
    Analyzes the amino acid patch and suggests mutations to relieve negative charge.
    """
    original_patch = {
        47: {"name": "Serine", "one_letter": "S"},
        48: {"name": "Glutamate", "one_letter": "E"},
        49: {"name": "Glutamate", "one_letter": "E"},
        50: {"name": "Aspartate", "one_letter": "D"}
    }

    mutant_amino_acid = {"name": "Alanine", "one_letter": "A"}

    print("### Plan for Site-Directed Mutagenesis ###\n")
    print("Objective: To relieve the negative charge of the S-E-E-D patch at positions 47-50.")
    print(f"Strategy: Replace all four residues with a neutral, non-phosphorylatable, and structurally conservative amino acid: {mutant_amino_acid['name']} ({mutant_amino_acid['one_letter']}).\n")

    original_sequence = []
    mutant_sequence = []

    print("--- Proposed Amino Acid Changes ---")
    for position in sorted(original_patch.keys()):
        original_aa = original_patch[position]
        print(f"Position {position}: Change original {original_aa['name']} ({original_aa['one_letter']}) to new {mutant_amino_acid['name']} ({mutant_amino_acid['one_letter']})")
        original_sequence.append(original_aa['one_letter'])
        mutant_sequence.append(mutant_amino_acid['one_letter'])
    
    print("\n--- Final Sequence Equation ---")
    print(f"Original sequence at 47-50: {' - '.join(original_sequence)}")
    print(f"Proposed mutant sequence at 47-50: {' - '.join(mutant_sequence)}")

if __name__ == '__main__':
    suggest_mutagenesis()
    # The final answer representing the new amino acid sequence from position 47 to 50
    final_answer = "AAAA"
    print(f"\n<<<{final_answer}>>>")
