def design_mutagenesis_experiment():
    """
    Analyzes the protein segment and proposes a site-directed mutagenesis
    experiment to eliminate its negative charge.
    """

    # Define the properties of the original amino acid patch
    original_patch = {
        47: {'code': 'S', 'name': 'Serine', 'property': 'Phosphorylatable, adding negative charge'},
        48: {'code': 'E', 'name': 'Glutamate', 'property': 'Acidic, negatively charged'},
        49: {'code': 'E', 'name': 'Glutamate', 'property': 'Acidic, negatively charged'},
        50: {'code': 'D', 'name': 'Aspartate', 'property': 'Acidic, negatively charged'}
    }

    # Define the proposed replacement
    replacement_aa = {'code': 'A', 'name': 'Alanine'}

    print("--- Analysis and Recommendation ---")
    print("The goal is to relieve the autoinhibitory negative charge of the patch from amino acid positions 47 to 50.")
    print("The original sequence is Ser-Glu-Glu-Asp (SEED). The negative charge comes from the acidic residues (E, D) and phosphorylation of Serine (S).\n")

    print(f"The best replacement amino acid for this purpose is {replacement_aa['name']} ({replacement_aa['code']}).")
    print("Rationale for mutating all four residues to Alanine:")
    print("1. Neutralization: Alanine is small, neutral, and chemically inert, effectively removing the negative charges from Glu and Asp.")
    print("2. Prevents Phosphorylation: Replacing Serine with Alanine removes the hydroxyl group, thereby preventing phosphorylation and the associated negative charge at position 47.")
    print("3. Gold Standard: This is a standard 'Alanine scanning' approach to systematically remove side-chain function without introducing major structural changes.\n")

    print("--- Final Mutation Equation ---")
    print("To create the mutant and test the hypothesis, the following mutations should be engineered.")
    print("Each mutation is listed in the standard 'OriginalAminoAcid-Position-NewAminoAcid' format:")

    # Construct and print the final "equation" by listing each individual mutation
    final_mutations = []
    for position in sorted(original_patch.keys()):
        original_code = original_patch[position]['code']
        new_code = replacement_aa['code']
        # The output format explicitly includes the position number for each part of the "equation"
        final_mutations.append(f"{original_code}{position}{new_code}")
    
    print(", ".join(final_mutations))

if __name__ == '__main__':
    design_mutagenesis_experiment()