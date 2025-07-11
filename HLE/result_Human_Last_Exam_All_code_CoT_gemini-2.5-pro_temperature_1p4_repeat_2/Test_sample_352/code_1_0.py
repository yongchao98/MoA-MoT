def design_mutagenesis():
    """
    Suggests and explains the best amino acid replacements to neutralize
    a negatively charged patch in protein x.
    """
    # Define the original amino acids and their positions
    original_patch = {
        47: {"aa": "Serine", "code": "S"},
        48: {"aa": "Glutamate", "code": "E"},
        49: {"aa": "Glutamate", "code": "E"},
        50: {"aa": "Aspartate", "code": "D"}
    }

    # Define the proposed replacement amino acid
    replacement_aa = {"aa": "Alanine", "code": "A"}

    # --- Print Explanation ---
    print("### Site-Directed Mutagenesis Plan ###\n")
    print("Objective: To relieve the inhibitory negative charge of the amino acid patch from position 47 to 50.\n")

    # Original patch details
    original_sequence = "-".join([info["code"] for pos, info in sorted(original_patch.items())])
    print(f"Original sequence at positions 47-50: {original_sequence}\n")
    print("Rationale for mutation:")
    print(f"- Position 47 ({original_patch[47]['aa']}): Is a phosphorylation site, which becomes highly negative when phosphorylated.")
    print(f"- Positions 48, 49 ({original_patch[48]['aa']}), and 50 ({original_patch[50]['aa']}): Are acidic and natively negatively charged.")
    print(f"\nProposed replacement: Replace all four residues with {replacement_aa['aa']} ({replacement_aa['code']}).\n")

    # Proposed mutation details
    print("Justification:")
    print(f"- {replacement_aa['aa']} is small, neutral, non-polar, and cannot be phosphorylated.")
    print("- This change will neutralize the entire patch, removing both the native negative charges and the potential for charge addition via phosphorylation.\n")

    print("--- Proposed Final Equation ---")
    mutant_notation = []
    for position, info in sorted(original_patch.items()):
        mutation = f"{info['code']}{position}{replacement_aa['code']}"
        mutant_notation.append(mutation)
        print(f"Mutation for position {position}: Replace {info['aa']} ({info['code']}) with {replacement_aa['aa']} ({replacement_aa['code']}). Notation: {mutation}")

    final_mutant_name = ",".join(mutant_notation)
    final_sequence = replacement_aa['code'] * 4
    print(f"\nThe final mutant is designated: {final_mutant_name}")
    print(f"The new sequence at positions 47-50 will be: {final_sequence}")


# Execute the function to print the plan
design_mutagenesis()
<<<AAAA>>>