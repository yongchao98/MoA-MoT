def design_mutagenesis_experiment():
    """
    Outlines a site-directed mutagenesis plan to relieve autoinhibition
    in protein x by neutralizing a negatively charged amino acid patch.
    """

    # Define the original amino acid patch responsible for the negative charge.
    original_patch = {
        47: {"code": "S", "name": "Serine"},
        48: {"code": "E", "name": "Glutamate"},
        49: {"code": "E", "name": "Glutamate"},
        50: {"code": "D", "name": "Aspartate"}
    }

    # Define the proposed replacement amino acids.
    # Alanine (A) is chosen to neutralize charge and remove the phosphorylation site
    # with minimal structural disruption.
    proposed_amino_acid = {"code": "A", "name": "Alanine"}

    print("--- Mutagenesis Plan to Relieve Protein X Autoinhibition ---")
    print("\nObjective: Neutralize the highly negative patch at amino acid positions 47-50.")
    print("Strategy: Replace the original residues with Alanine, a small and neutral amino acid.\n")

    print("1. Original Amino Acid Patch (S-E-E-D):")
    for position, aa_info in sorted(original_patch.items()):
        print(f"   Position {position}: {aa_info['name']} ({aa_info['code']}) - This residue contributes to the negative charge (or can be phosphorylated).")

    print("\n2. Proposed Mutations (A-A-A-A):")
    final_mutation_str_parts = []
    for position, aa_info in sorted(original_patch.items()):
        mutation_notation = f"{aa_info['code']}{position}{proposed_amino_acid['code']}"
        final_mutation_str_parts.append(mutation_notation)
        print(f"   At Position {position}: Change {aa_info['name']} ({aa_info['code']}) to {proposed_amino_acid['name']} ({proposed_amino_acid['code']}).")

    final_mutation_equation = ", ".join(final_mutation_str_parts)
    print("\n--------------------------------------------------------------")
    print("Final Recommended Set of Mutations (Standard Notation):")
    # This loop outputs each number in the final mutation notation as requested
    print_parts = []
    for part in final_mutation_str_parts:
        original_aa = part[0]
        position_num = ''.join(filter(str.isdigit, part))
        new_aa = part[-1]
        print_parts.append(f"{original_aa}{position_num}{new_aa}")
    
    # We will print the final equation/string showing all mutations together
    print(f"The final proposed mutation set is: {final_mutation_equation}")
    print("--------------------------------------------------------------")
    print("\nRationale: This S47A, E48A, E49A, D50A mutant will lack the negative charges and")
    print("the phosphorylation site, directly testing the hypothesis that this specific")
    print("patch is responsible for the autoinhibitory activity of the IDR.")


# Execute the function to display the plan.
if __name__ == "__main__":
    design_mutagenesis_experiment()

<<<The best replacement for the S-E-E-D patch at positions 47-50 is A-A-A-A, corresponding to the mutations S47A, E48A, E49A, and D50A.>>>