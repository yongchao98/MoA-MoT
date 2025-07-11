def suggest_mutagenesis():
    """
    Suggests and explains the best amino acid replacements for a site-directed
    mutagenesis experiment to relieve an inhibitory negatively charged patch.
    """
    
    # Define the original amino acid patch
    original_patch = {
        47: {'name': 'Serine', 'code': 'S'},
        48: {'name': 'Glutamate', 'code': 'E'},
        49: {'name': 'Glutamate', 'code': 'E'},
        50: {'name': 'Aspartate', 'code': 'D'}
    }

    # Define the proposed mutant patch
    # The best replacement is Alanine (A) for all positions.
    mutant_patch = {
        47: {'name': 'Alanine', 'code': 'A'},
        48: {'name': 'Alanine', 'code': 'A'},
        49: {'name': 'Alanine', 'code': 'A'},
        50: {'name': 'Alanine', 'code': 'A'}
    }

    print("--- Site-Directed Mutagenesis Plan ---")
    print("\nObjective: To relieve the autoinhibitory effect of the negatively charged S-E-E-D patch (positions 47-50).")

    print("\nRationale for Alanine Replacement:")
    print("1. Neutralization of Charge: Alanine is a neutral amino acid. Replacing the negatively charged Glutamate (E) and Aspartate (D) will eliminate the patch's negative charge.")
    print("2. Prevention of Phosphorylation: Replacing Serine (S) with Alanine (A) removes the hydroxyl group required for phosphorylation, thus preventing the introduction of a strong negative charge at position 47.")
    print("3. Minimal Structural Impact: Alanine is a small amino acid, which minimizes the risk of disrupting the protein's overall structure.")

    print("\n--- Proposed Mutation Details ---")
    original_seq_str = "".join([info['code'] for pos, info in sorted(original_patch.items())])
    mutant_seq_str = "".join([info['code'] for pos, info in sorted(mutant_patch.items())])
    
    print(f"Original sequence at 47-50: {original_seq_str}")
    print(f"Proposed mutant sequence at 47-50: {mutant_seq_str}\n")
    
    print("Specific Amino Acid Changes:")
    for pos in sorted(original_patch.keys()):
        orig_aa = original_patch[pos]
        mut_aa = mutant_patch[pos]
        print(f"Position {pos}: Replace {orig_aa['name']} ({orig_aa['code']}) with {mut_aa['name']} ({mut_aa['code']})")

suggest_mutagenesis()

print("\n<<<AAAA>>>")