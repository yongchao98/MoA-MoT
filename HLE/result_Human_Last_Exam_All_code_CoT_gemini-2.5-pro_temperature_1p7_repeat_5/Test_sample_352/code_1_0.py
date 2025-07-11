def suggest_mutagenesis():
    """
    Analyzes an amino acid patch and suggests mutations to remove its negative charge.
    """
    original_patch = {
        47: "S (Serine)",
        48: "E (Glutamate)",
        49: "E (Glutamate)",
        50: "D (Aspartate)"
    }

    replacement_amino_acid = "A (Alanine)"

    print("--- Experimental Design for Site-Directed Mutagenesis ---")
    print("\nObjective: Relieve the negative charge from the S-E-E-D patch at positions 47-50.\n")
    print("Rationale:")
    print("The goal is to test the hypothesis that the patch's negative charge causes autoinhibition.")
    print("To do this, we replace the original amino acids with a small, neutral, and non-polar amino acid.")
    print(f"The ideal choice is {replacement_amino_acid}.\n")
    print("1. Replacing Serine (S47) with Alanine (A) prevents phosphorylation, thereby eliminating a key source of negative charge.")
    print("2. Replacing negatively charged Glutamate (E48, E49) and Aspartate (D50) with neutral Alanine removes the acidic side chains.\n")
    
    print("--- Recommended Mutations ---")
    original_sequence_str = ""
    mutant_sequence_str = ""
    for position, aa in sorted(original_patch.items()):
        original_aa_code = aa.split(" ")[0]
        replacement_aa_code = replacement_amino_acid.split(" ")[0]
        
        # Build the string for the equation
        if original_sequence_str:
            original_sequence_str += "-"
            mutant_sequence_str += "-"
        original_sequence_str += f"{original_aa_code}{position}"
        mutant_sequence_str += f"{replacement_aa_code}{position}"
        
        # Print individual mutation
        print(f"At position {position}: Mutate {aa} to {replacement_amino_acid}.  (Notation: {original_aa_code}{position}{replacement_aa_code})")

    print("\n--- Final Equation ---")
    print("The final transformation replaces the entire inhibitory patch:")
    print(f"{original_sequence_str}  ->  {mutant_sequence_str}")


if __name__ == '__main__':
    suggest_mutagenesis()
    # The following line is a placeholder for the final answer as per instructions.
    # The code above provides the detailed reasoning and recommended mutations.
    # The best replacement is to mutate all four residues to Alanine (A).
    print("\n<<<S47A, E48A, E49A, D50A>>>")