def design_mutagenesis():
    """
    Designs and explains the site-directed mutagenesis experiment to neutralize
    a negatively charged patch in protein x.
    """

    # Define the original amino acids and their positions
    original_patch = {
        47: {"one_letter_code": "S", "full_name": "Serine"},
        48: {"one_letter_code": "E", "full_name": "Glutamate"},
        49: {"one_letter_code": "E", "full_name": "Glutamate"},
        50: {"one_letter_code": "D", "full_name": "Aspartate"}
    }

    # Define the replacement amino acid
    replacement_aa = {"one_letter_code": "A", "full_name": "Alanine"}

    # --- Explanation ---
    print("Experimental Design: Site-Directed Mutagenesis of Protein X")
    print("Objective: To relieve the autoinhibitory negative charge of the patch at positions 47-50.\n")
    print("Rationale:")
    print("The original patch S-E-E-D is highly negatively charged due to the acidic residues (E, D) and the phosphorylation of Serine (S).")
    print("To eliminate this negative charge, each residue will be replaced with Alanine (A). Alanine is a small, neutral, and non-phosphorylatable amino acid, making it the ideal choice to test the function of the charged side chains.\n")

    # --- Proposed Mutation ---
    print("Proposed Mutations:")
    original_sequence = []
    mutated_sequence = []
    
    for position, aa_info in sorted(original_patch.items()):
        original_sequence.append(aa_info['one_letter_code'])
        mutated_sequence.append(replacement_aa['one_letter_code'])
        print(f"Position {position}: Replace original {aa_info['one_letter_code']} ({aa_info['full_name']}) with {replacement_aa['one_letter_code']} ({replacement_aa['full_name']})")

    original_str = "".join(original_sequence)
    mutated_str = "".join(mutated_sequence)
    
    print(f"\nSummary: The original amino acid sequence '{original_str}' at positions 47-50 should be mutated to '{mutated_str}'.")

if __name__ == '__main__':
    design_mutagenesis()
    # The final answer in the required format
    final_answer = "The replacement of the amino acid sequence Serine-Glutamate-Glutamate-Aspartate (SEED) at positions 47-50 with Alanine-Alanine-Alanine-Alanine (AAAA)"
    print(f"\n<<<{final_answer}>>>")
