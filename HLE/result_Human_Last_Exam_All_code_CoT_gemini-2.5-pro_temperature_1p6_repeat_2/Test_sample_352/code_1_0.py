def suggest_mutagenesis():
    """
    Suggests and explains the best amino acid replacements for a site-directed
    mutagenesis experiment to neutralize a negatively charged protein patch.
    """
    original_patch = {
        47: {"aa": "S", "name": "Serine"},
        48: {"aa": "E", "name": "Glutamate"},
        49: {"aa": "E", "name": "Glutamate"},
        50: {"aa": "D", "name": "Aspartate"}
    }

    replacement_aa = {"aa": "A", "name": "Alanine"}

    print("### Site-Directed Mutagenesis Plan ###\n")
    print(f"Goal: To relieve the autoinhibitory effect of the negatively charged patch at amino acid positions 47-50.\n")
    print(f"Original Sequence: ...-{original_patch[47]['aa']}{original_patch[48]['aa']}{original_patch[49]['aa']}{original_patch[50]['aa']}-... (at positions 47-50)\n")

    print("--- Analysis and Recommendation ---\n")

    # Position 47
    pos = 47
    orig = original_patch[pos]
    print(f"Position {pos}: Original is {orig['name']} ({orig['aa']})")
    print(f"  - Problem: This is a phosphorylation site. Phosphorylation adds a strong negative charge.")
    print(f"  - Solution: Replace with {replacement_aa['name']} ({replacement_aa['aa']}).")
    print(f"  - Justification: Alanine is neutral and cannot be phosphorylated, thus eliminating the potential negative charge at this position.")
    print(f"  - Proposed Mutation: {pos}{orig['aa']} -> {pos}{replacement_aa['aa']}\n")

    # Positions 48, 49, 50
    for pos in [48, 49, 50]:
        orig = original_patch[pos]
        print(f"Position {pos}: Original is {orig['name']} ({orig['aa']})")
        print(f"  - Problem: This is an acidic amino acid, carrying a negative charge at physiological pH.")
        print(f"  - Solution: Replace with {replacement_aa['name']} ({replacement_aa['aa']}).")
        print(f"  - Justification: Alanine is small, neutral, and non-polar, which effectively removes the negative charge with minimal structural disruption.")
        print(f"  - Proposed Mutation: {pos}{orig['aa']} -> {pos}{replacement_aa['aa']}\n")

    final_original_sequence = "".join([original_patch[p]['aa'] for p in sorted(original_patch.keys())])
    final_mutated_sequence = replacement_aa['aa'] * 4

    print("--- Summary ---")
    print(f"The recommended experiment is to mutate the original sequence '{final_original_sequence}' at positions 47-50")
    print(f"to the new sequence: '{final_mutated_sequence}'.")
    print("This tetra-Alanine (AAAA) mutant will create a neutral, non-phosphorylatable patch,")
    print("allowing for a clear test of the hypothesis that the negative charge of this region is autoinhibitory.")

if __name__ == '__main__':
    suggest_mutagenesis()