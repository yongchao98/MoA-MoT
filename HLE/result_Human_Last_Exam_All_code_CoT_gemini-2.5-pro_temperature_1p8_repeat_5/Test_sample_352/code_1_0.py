def suggest_mutagenesis():
    """
    Suggests amino acid replacements for a site-directed mutagenesis experiment
    to relieve a negatively charged patch in protein x.
    """

    # Define the original amino acid patch
    positions = [47, 48, 49, 50]
    original_aa = ['Serine (S)', 'Glutamate (E)', 'Glutamate (E)', 'Aspartate (D)']
    original_charges = ['Can be negatively charged (via phosphorylation)', 'Negative', 'Negative', 'Negative']

    # Define the proposed replacement amino acid
    replacement_aa = 'Alanine (A)'
    replacement_charge = 'Neutral'

    print("--- Site-Directed Mutagenesis Plan ---")
    print("\nObjective: To relieve the autoinhibitory negative charge of the patch from positions 47-50.\n")

    print("Original Patch Analysis:")
    print("------------------------")
    for i in range(len(positions)):
        print(f"Position {positions[i]}: {original_aa[i]:<14} | Charge: {original_charges[i]}")
    print("\nThis patch (S-E-E-D) is highly negatively charged, contributing to autoinhibition.")

    print("\nProposed Mutation Strategy:")
    print("--------------------------")
    print("Replace all four residues with Alanine (A). Alanine is small, neutral, and minimally disruptive.")
    print("This 'Alanine Scan' will neutralize the patch and prevent phosphorylation at position 47.\n")

    print("Proposed Mutant Patch:")
    print("-----------------------")
    for i in range(len(positions)):
        print(f"Position {positions[i]}: Original {original_aa[i].split(' ')[0]} -> Proposed {replacement_aa.split(' ')[0]}")

    print("\nFinal Recommended Mutant Sequence at positions 47-50: A-A-A-A")
    print("(Mutations: S47A, E48A, E49A, D50A)")

if __name__ == '__main__':
    suggest_mutagenesis()