def design_mutagenesis():
    """
    Prints the plan for a site-directed mutagenesis experiment to neutralize
    a negatively charged patch in a protein.
    """

    # Define the original and proposed mutant amino acids for the patch
    patch = {
        47: {'original': 'Serine (S)', 'mutant': 'Alanine (A)', 'rationale': 'Remove phosphorylation site to prevent negative charge addition.'},
        48: {'original': 'Glutamate (E)', 'mutant': 'Alanine (A)', 'rationale': 'Replace negatively charged residue with a neutral one.'},
        49: {'original': 'Glutamate (E)', 'mutant': 'Alanine (A)', 'rationale': 'Replace negatively charged residue with a neutral one.'},
        50: {'original': 'Aspartate (D)', 'mutant': 'Alanine (A)', 'rationale': 'Replace negatively charged residue with a neutral one.'}
    }

    print("### Mutagenesis Plan to Neutralize the Inhibitory Patch ###\n")
    print("The optimal strategy is to replace each amino acid in the 47-50 patch with Alanine (A).")
    print("This approach, known as Alanine scanning, effectively eliminates the negative charge and the potential for phosphorylation with minimal structural disruption.\n")

    print("--- Detailed Mutations ---")
    original_sequence = []
    mutant_sequence = []

    for position in sorted(patch.keys()):
        info = patch[position]
        print(f"Position {position}:")
        print(f"  - Replace Original: {info['original']}")
        print(f"  - With Mutant:      {info['mutant']}")
        print(f"  - Rationale:        {info['rationale']}\n")
        original_sequence.append(info['original'].split(' ')[1].replace('(', '').replace(')', ''))
        mutant_sequence.append(info['mutant'].split(' ')[1].replace('(', '').replace(')', ''))

    print("--- Summary of Change ---")
    print(f"Original Sequence Patch (47-50): {original_sequence[0]}-{original_sequence[1]}-{original_sequence[2]}-{original_sequence[3]}")
    print(f"Proposed Mutant Sequence Patch:    {mutant_sequence[0]}-{mutant_sequence[1]}-{mutant_sequence[2]}-{mutant_sequence[3]}")
    print("\nThis experiment creates the S47A/E48A/E49A/D50A quadruple mutant.")

if __name__ == '__main__':
    design_mutagenesis()
<<<AAAA>>>