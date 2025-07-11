def suggest_mutagenesis():
    """
    Suggests and explains the best amino acid replacements for a site-directed
    mutagenesis experiment to relieve a negatively charged patch in protein x.
    """
    # Define the original amino acids and their positions
    original_patch = {
        47: 'Serine (S)',
        48: 'Glutamate (E)',
        49: 'Glutamate (E)',
        50: 'Aspartate (D)'
    }

    # The chosen replacement amino acid
    replacement_aa = 'Alanine (A)'

    print("### Site-Directed Mutagenesis Plan ###\n")
    print("Goal: To neutralize the highly negative charge of the patch from amino acid 47 to 50.\n")
    print(f"Original Patch: {original_patch[47]}-{original_patch[48]}-{original_patch[49]}-{original_patch[50]}\n")
    print(f"Proposed Replacement: The best strategy is to replace all four residues with {replacement_aa}.\n")

    print("Rationale:")
    print(f"1. Neutralization: Alanine is a small, neutral, and non-polar amino acid. Replacing the negatively charged Glutamate (E) and Aspartate (D) with Alanine will remove the negative charges.")
    print(f"2. Preventing Phosphorylation: Serine (S) at position 47 is a phosphorylation site, which introduces a strong negative charge. Replacing Serine with Alanine (S47A) removes the hydroxyl group required for phosphorylation, thus preventing this modification.")
    print(f"3. Minimal Disruption: This 'Alanine scan' is a standard method to test the function of specific residues while minimizing the risk of introducing major structural changes or new, unintended chemical interactions.\n")

    print("--- Recommended Mutations ---")
    print("The original sequence S-E-E-D at positions 47-50 should be mutated to A-A-A-A.")
    print("The specific mutations to engineer are:")

    # Print the final equation/mutations
    mutations = []
    for position, aa_name in original_patch.items():
        original_code = aa_name.split('(')[1][0]
        mutation_string = f"Original: {position} {aa_name} -> New: {position} {replacement_aa}"
        mutations.append(f"{original_code}{position}A")
        print(mutation_string)

    print(f"\nFinal mutation set in standard notation: {', '.join(mutations)}")

suggest_mutagenesis()
<<<S47A, E48A, E49A, D50A>>>