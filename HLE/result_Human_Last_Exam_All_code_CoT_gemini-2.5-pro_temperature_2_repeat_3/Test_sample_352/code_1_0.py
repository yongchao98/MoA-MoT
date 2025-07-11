def suggest_mutations():
    """
    Suggests and explains amino acid mutations to relieve the negative charge
    of a specific protein patch.
    """

    # Original amino acid sequence at the specified positions
    original_sequence = {
        47: {"code": "S", "name": "Serine"},
        48: {"code": "E", "name": "Glutamate"},
        49: {"code": "E", "name": "Glutamate"},
        50: {"code": "D", "name": "Aspartate"}
    }

    # Proposed replacements to create a neutral patch
    # S -> A (Alanine): Prevents phosphorylation, neutral.
    # E -> Q (Glutamine): Neutral amide of Glutamate.
    # D -> N (Asparagine): Neutral amide of Aspartate.
    proposed_mutations = {
        47: {"code": "A", "name": "Alanine"},
        48: {"code": "Q", "name": "Glutamine"},
        49: {"code": "Q", "name": "Glutamine"},
        50: {"code": "N", "name": "Asparagine"}
    }

    print("Experiment Design: Site-Directed Mutagenesis of Protein X Patch (47-50)")
    print("-" * 70)
    print("Goal: To relieve the inhibitory negative charge of the S-E-E-D sequence.\n")

    print("The proposed mutations are designed to neutralize the region while conserving structural similarity:\n")

    final_sequence_code = ""
    for position in sorted(original_sequence.keys()):
        original_aa = original_sequence[position]
        mutant_aa = proposed_mutations[position]
        final_sequence_code += mutant_aa["code"]

        print(
            f"Position {position}: Replace original {original_aa['code']} ({original_aa['name']}) "
            f"with proposed {mutant_aa['code']} ({mutant_aa['name']})."
        )

    original_code_string = "".join(original_sequence[pos]['code'] for pos in sorted(original_sequence.keys()))

    print("\n--- Summary ---")
    print(f"The native inhibitory sequence at positions 47-50 is: {original_code_string}")
    print(f"The proposed neutral mutant sequence is: {final_sequence_code}")
    print(
        "\nThis change eliminates the phosphorylation site and replaces negatively charged residues with their neutral counterparts, "
        "providing the best test for the autoinhibition hypothesis."
    )
    print(f"<<<{final_sequence_code}>>>")


if __name__ == "__main__":
    suggest_mutations()