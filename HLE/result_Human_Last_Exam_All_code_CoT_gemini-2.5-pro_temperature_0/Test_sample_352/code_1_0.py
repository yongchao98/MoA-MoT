def suggest_mutagenesis():
    """
    Suggests and explains a site-directed mutagenesis experiment
    to neutralize a negatively charged protein patch.
    """
    # Define the original amino acids and their positions
    original_patch = {
        47: {'code': 'S', 'name': 'Serine'},
        48: {'code': 'E', 'name': 'Glutamate'},
        49: {'code': 'E', 'name': 'Glutamate'},
        50: {'code': 'D', 'name': 'Aspartate'}
    }

    # Define the proposed replacement amino acid
    # Alanine (A) is chosen because it is small, neutral, non-polar,
    # and cannot be phosphorylated. This effectively removes the negative
    # charge of the patch with minimal structural disruption.
    replacement_aa = {'code': 'A', 'name': 'Alanine'}

    # Print the rationale and the proposed mutation plan
    print("Experimental Design: Site-Directed Mutagenesis of Protein X")
    print("="*60)
    print("Objective: To eliminate the autoinhibitory negative charge of the patch at positions 47-50.")
    print("\nRationale: The original S-E-E-D sequence is highly negative due to the acidic residues (E, D) and a potential phosphorylation site (S). Replacing these with a small, neutral amino acid like Alanine will test the hypothesis that this negative charge is responsible for autoinhibition.")
    print("\nProposed Mutation Plan:")
    print("-" * 60)
    print(f"{'Position':<10} | {'Original Amino Acid':<25} | {'Proposed Mutant':<20}")
    print("-" * 60)

    original_sequence = []
    mutated_sequence = []

    for position, aa_info in original_patch.items():
        original_name = f"{aa_info['name']} ({aa_info['code']})"
        mutant_name = f"{replacement_aa['name']} ({replacement_aa['code']})"
        print(f"{position:<10} | {original_name:<25} | {mutant_name:<20}")
        original_sequence.append(aa_info['code'])
        mutated_sequence.append(replacement_aa['code'])

    print("-" * 60)
    print("\nSummary of Change:")
    print(f"The original sequence from position 47 to 50: {'-'.join(original_sequence)}")
    print(f"Will be mutated to the new sequence:           {'-'.join(mutated_sequence)}")
    print("="*60)

# Execute the function to print the plan
suggest_mutagenesis()