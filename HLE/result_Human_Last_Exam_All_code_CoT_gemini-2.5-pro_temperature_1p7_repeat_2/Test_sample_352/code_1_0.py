def suggest_mutagenesis():
    """
    Analyzes the protein sequence patch and suggests the best replacements
    for a site-directed mutagenesis experiment based on the provided context.
    """
    # Define the original amino acid sequence patch
    original_sequence = {
        47: {"code": "S", "name": "Serine"},
        48: {"code": "E", "name": "Glutamate"},
        49: {"code": "E", "name": "Glutamate"},
        50: {"code": "D", "name": "Aspartate"}
    }

    # Define the rationale for the negative charge of the patch
    charge_properties = {
        47: "Serine (S) is a phosphorylation site, becoming highly negative when phosphorylated.",
        48: "Glutamate (E) is an acidic amino acid, negatively charged at physiological pH.",
        49: "Glutamate (E) is an acidic amino acid, negatively charged at physiological pH.",
        50: "Aspartate (D) is an acidic amino acid, negatively charged at physiological pH."
    }

    # Define the proposed mutations and the rationale
    proposed_mutations = {
        47: {"code": "A", "name": "Alanine", "rationale": "Cannot be phosphorylated, removing the potential for a negative charge."},
        48: {"code": "Q", "name": "Glutamine", "rationale": "Structurally similar to Glutamate but is neutral, directly testing the role of the charge."},
        49: {"code": "Q", "name": "Glutamine", "rationale": "Structurally similar to Glutamate but is neutral, directly testing the role of the charge."},
        50: {"code": "N", "name": "Asparagine", "rationale": "Structurally similar to Aspartate but is neutral, directly testing the role of the charge."}
    }

    print("--- Analysis of Mutagenesis Target ---")
    original_seq_str = "-".join([d["code"] for d in original_sequence.values()])
    print(f"Original sequence at positions 47-50: {original_seq_str}\n")
    print("This patch is inhibitory due to a high negative charge from its amino acids:")
    for pos, reason in charge_properties.items():
        print(f"  - Position {pos} ({original_sequence[pos]['code']}): {reason}")

    print("\n--- Proposed Mutations to Neutralize the Patch ---")
    mutant_seq_str = "-".join([d["code"] for d in proposed_mutations.values()])
    print("To best test the hypothesis that the negative charge causes inhibition,")
    print("the following conservative, charge-neutralizing mutations are recommended:\n")

    for pos, mutation_info in proposed_mutations.items():
        original_aa = f"{original_sequence[pos]['name']} ({original_sequence[pos]['code']})"
        mutant_aa = f"{mutation_info['name']} ({mutation_info['code']})"
        print(f"Position {pos}: Mutate {original_aa} to {mutant_aa}")
        print(f"  Rationale: {mutation_info['rationale']}")

    print("\n--- Final Recommended Mutant Sequence ---")
    print(f"The final recommended sequence is: {mutant_seq_str} (S47A, E48Q, E49Q, D50N)")

if __name__ == '__main__':
    suggest_mutagenesis()