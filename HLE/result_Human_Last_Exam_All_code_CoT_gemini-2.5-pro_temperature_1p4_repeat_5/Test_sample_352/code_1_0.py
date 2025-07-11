import sys

def design_mutagenesis():
    """
    Designs and explains a site-directed mutagenesis experiment to relieve
    an inhibitory, negatively charged amino acid patch.
    """
    # --- Define Original and Proposed Mutant Sequences ---
    original_patch = {
        47: {'short': 'S', 'long': 'Serine'},
        48: {'short': 'E', 'long': 'Glutamate'},
        49: {'short': 'E', 'long': 'Glutamate'},
        50: {'short': 'D', 'long': 'Aspartate'}
    }
    replacement_amino_acid = {'short': 'A', 'long': 'Alanine'}

    # --- Print Header ---
    print("Plan for Site-Directed Mutagenesis of Protein X")
    print("=" * 50)
    print("\nObjective: To eliminate the autoinhibitory effect of the negatively charged\npatch at amino acid positions 47-50.\n")

    # --- Print Analysis of the Original Patch ---
    print("1. Analysis of the Original Patch (S-E-E-D):")
    original_sequence_str = '-'.join([f"{aa_info['short']}{pos}" for pos, aa_info in original_patch.items()])
    print(f"   The original sequence is: {original_sequence_str}")
    print(f"   - {original_patch[47]['long']} (S47): A potential phosphorylation site. When phosphorylated, it adds a strong negative charge.")
    print(f"   - {original_patch[48]['long']} (E48) & {original_patch[49]['long']} (E49): Acidic and negatively charged.")
    print(f"   - {original_patch[50]['long']} (D50): Acidic and negatively charged.\n")

    # --- Print Proposed Solution ---
    print(f"2. Proposed Replacement Amino Acid: {replacement_amino_acid['long']} ({replacement_amino_acid['short']})")
    print("   Rationale for choosing Alanine:")
    print("     a) Eliminates Negative Charge: Alanine is a neutral amino acid, which will remove the negative charges from E48, E49, and D50.")
    print("     b) Prevents Phosphorylation: Replacing Serine with Alanine removes the hydroxyl group required for phosphorylation, thereby preventing any charge addition at position 47.")
    print("     c) Minimally Disruptive: Alanine is small and biochemically inert, making it the standard choice for testing the functional role of amino acid side chains without introducing major structural changes.\n")

    # --- Print the Final Proposed Mutation "Equation" ---
    print("3. Final Proposed Mutation:")
    mutant_sequence_str = '-'.join([f"{replacement_amino_acid['short']}{pos}" for pos in original_patch.keys()])
    
    # Building the final equation string showing each individual mutation
    mutation_list = []
    for pos, aa_info in original_patch.items():
        mutation_list.append(f"{aa_info['short']}{pos}{replacement_amino_acid['short']}")
    
    print("   The proposed series of mutations is to change the 'SEED' patch to 'AAAA'.")
    print("\n   Mutation Equation:")
    print(f"   Original Patch: {original_patch[47]['short']}47 - {original_patch[48]['short']}48 - {original_patch[49]['short']}49 - {original_patch[50]['short']}50")
    print(f"   Proposed Mutant:  {replacement_amino_acid['short']}47 - {replacement_amino_acid['short']}48 - {replacement_amino_acid['short']}49 - {replacement_amino_acid['short']}50")

    # --- Prepare the final answer for the <<<>>> format ---
    # Standard mutation notation is e.g., S47A
    final_answer = ", ".join(mutation_list)
    # This line will be captured to generate the final answer block.
    # It won't be visible in the script's regular output.
    sys.stdout.write(f"\n<<<{final_answer}>>>")


if __name__ == '__main__':
    design_mutagenesis()