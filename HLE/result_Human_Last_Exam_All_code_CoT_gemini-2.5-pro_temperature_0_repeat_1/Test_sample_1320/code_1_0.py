def find_helix_type():
    """
    Determines the likely helix type for a hybrid peptide based on the average
    backbone length of its constituent residues.
    """
    # Step 1 & 2: Define the number of backbone atoms for each residue type.
    # Alanine is an alpha-amino acid.
    alpha_atoms = 3  # Backbone path: N-C(alpha)-C'
    # The other residue is an epsilon-amino acid.
    epsilon_atoms = 7 # Backbone path: N-C(epsilon)-...-C(alpha)-C'

    print("Step 1: Identify the number of backbone atoms for each residue type.")
    print(f"An alpha-amino acid (alanine) has {alpha_atoms} backbone atoms.")
    print(f"An epsilon-amino acid has {epsilon_atoms} backbone atoms.")
    print("-" * 20)

    # Step 3: Calculate the average number of backbone atoms for the alternating sequence.
    avg_atoms = (alpha_atoms + epsilon_atoms) / 2

    print("Step 2: Calculate the average backbone atoms for the alternating (alpha-epsilon) sequence.")
    print(f"Average = ({alpha_atoms} + {epsilon_atoms}) / 2 = {int(avg_atoms)}")
    print("-" * 20)

    # Step 4 & 5: Use analogy to map the average length to a known helix type.
    # This mapping is based on established foldamer literature.
    helix_type_map = {
        4: "12/14",  # Beta-peptides
        5: "14/16",  # Gamma-peptides
        6: "16/18",  # Delta-peptides
        7: "18/20"   # Epsilon-peptides
    }

    print("Step 3: Apply the principle of analogy.")
    print(f"An average of {int(avg_atoms)} backbone atoms is equivalent to a gamma-peptide.")
    
    if int(avg_atoms) in helix_type_map:
        result_helix = helix_type_map[int(avg_atoms)]
        # The problem asks to output each number in the final equation.
        num1, num2 = result_helix.split('/')
        
        print(f"Homopolymers of gamma-peptides are known to form a characteristic {result_helix}-helix.")
        print("Therefore, the alpha/epsilon-hybrid foldamer is most likely to form this type of helix.")
        print("-" * 20)
        print("The final equation representing the helix type is composed of the following numbers:")
        print(num1)
        print(num2)
    else:
        print("Could not determine the helix type for this average backbone length.")

find_helix_type()