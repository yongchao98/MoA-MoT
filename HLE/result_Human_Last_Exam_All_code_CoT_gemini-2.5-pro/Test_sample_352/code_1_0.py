def design_mutagenesis_experiment():
    """
    Designs and explains a site-directed mutagenesis experiment to neutralize
    a negatively charged patch in a protein.
    """

    # Define the original amino acids and their positions in the protein.
    original_patch = {
        47: {"name": "Serine", "code": "S"},
        48: {"name": "Glutamate", "code": "E"},
        49: {"name": "Glutamate", "code": "E"},
        50: {"name": "Aspartate", "code": "D"}
    }

    # The best replacement amino acid is Alanine (A).
    replacement_aa = {"name": "Alanine", "code": "A"}

    print("### Site-Directed Mutagenesis Experiment Design ###")
    print("\nObjective: To eliminate the inhibitory negative charge from the patch at amino acid positions 47-50.")
    
    print("\nRationale:")
    print("The chosen replacement for all positions is Alanine (A). This is the standard for functional testing because:")
    print("1. It is electrically neutral, which will remove the negative charges from Glutamate (E) and Aspartate (D).")
    print("2. It cannot be phosphorylated, which removes the potential negative charge from Serine (S) at position 47.")
    print("3. Its small size minimizes the risk of disrupting the protein's structure.")

    print("\nProposed Mutation Details:")
    print("-" * 30)
    print(f"{'Position':<10} | {'Original AA':<15} | {'Proposed AA':<15}")
    print("-" * 30)

    for position, aa_info in original_patch.items():
        original_str = f"{aa_info['name']} ({aa_info['code']})"
        replacement_str = f"{replacement_aa['name']} ({replacement_aa['code']})"
        print(f"{position:<10} | {original_str:<15} | {replacement_str:<15}")
    
    print("-" * 30)
    original_sequence = "".join([info['code'] for pos, info in sorted(original_patch.items())])
    new_sequence = replacement_aa['code'] * 4
    print(f"\nSummary: Mutate the original sequence '{original_sequence}' (positions 47-50) to '{new_sequence}'.")


if __name__ == '__main__':
    design_mutagenesis_experiment()