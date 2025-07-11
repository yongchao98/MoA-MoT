def design_mutagenesis_experiment():
    """
    Designs and explains the site-directed mutagenesis experiment to relieve
    the autoinhibitory effect of a negatively charged amino acid patch.
    """

    # Define the original amino acids and their positions
    original_patch = {
        47: {"name": "Serine", "code": "S"},
        48: {"name": "Glutamate", "code": "E"},
        49: {"name": "Glutamate", "code": "E"},
        50: {"name": "Aspartate", "code": "D"}
    }

    # The chosen replacement amino acid
    replacement_aa = {"name": "Alanine", "code": "A"}

    print("--- Site-Directed Mutagenesis Plan ---")
    print(f"Objective: To neutralize the negatively charged patch at amino acid positions 47-50.\n")
    print("Strategy: Replace the original residues with a small, neutral amino acid to eliminate charge and prevent phosphorylation without disrupting the protein structure. The ideal candidate is Alanine (A).\n")

    print("Recommended Mutations:")
    mutated_sequence_code = ""
    for position, aa_info in original_patch.items():
        print(f"Position {position}: Replace {aa_info['name']} ({aa_info['code']}) with {replacement_aa['name']} ({replacement_aa['code']})")
        mutated_sequence_code += replacement_aa['code']

    print("\nSummary of the mutation:")
    original_sequence_str = "".join([info['code'] for info in original_patch.values()])
    print(f"Original sequence segment (47-50): {original_sequence_str}")
    print(f"Proposed mutated sequence segment (47-50): {mutated_sequence_code}")
    print("\nThis A-A-A-A mutation effectively removes the negative charges from Glutamate and Aspartate and eliminates the phosphorylation potential of Serine, directly testing the hypothesis of charge-based autoinhibition.")

    # Final answer in the specified format
    final_answer = "AAAA"
    print(f"\n<<<The proposed replacement sequence for S-E-E-D at positions 47-50 is {final_answer}>>>")


if __name__ == '__main__':
    design_mutagenesis_experiment()