def suggest_mutations():
    """
    Suggests and prints the best amino acid replacements to neutralize a negatively charged patch.
    """
    
    # Original patch: Amino Acid and Position
    original_patch = {
        47: 'S',  # Serine (can be phosphorylated to become negative)
        48: 'E',  # Glutamate (negative)
        49: 'E',  # Glutamate (negative)
        50: 'D'   # Aspartate (negative)
    }

    # Rationale for replacement:
    # S -> A (Alanine): Prevents phosphorylation, neutral charge.
    # E -> Q (Glutamine): Neutralizes charge, structurally similar.
    # D -> N (Asparagine): Neutralizes charge, structurally similar.
    
    # Proposed replacement patch
    replacement_patch = {
        47: 'A',
        48: 'Q',
        49: 'Q',
        50: 'N'
    }

    print("To neutralize the negatively charged patch from positions 47 to 50, the following site-directed mutations are recommended:")
    print("-" * 50)
    
    # Print each recommended mutation
    for position in sorted(original_patch.keys()):
        original_aa = original_patch[position]
        replacement_aa = replacement_patch[position]
        # Printing each part of the mutation "equation" as requested
        print(f"Position {position}: Replace original '{original_aa}' with new '{replacement_aa}'. Mutation: {original_aa}{position}{replacement_aa}")
        
    print("-" * 50)
    final_sequence = "".join(replacement_patch[pos] for pos in sorted(replacement_patch.keys()))
    print(f"The original S-E-E-D motif should be replaced with: {final_sequence}")

if __name__ == "__main__":
    suggest_mutations()
