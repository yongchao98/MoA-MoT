def calculate_basis_set_primitives():
    """
    Calculates and explains the number of primitive Gaussians for a Carbon atom
    in the 6-311G** basis set.
    """
    atom_name = "Carbon (C)"
    print(f"Calculating the number of primitive Gaussians for {atom_name} in a 6-311G** basis set.")
    print("---------------------------------------------------------------------------------")
    print("The 6-311G** basis set notation is broken down as follows:")
    print(" - '6-':       Core orbitals are described by 1 function built from 6 primitives.")
    print(" - '311':      Valence orbitals are split into 3 functions, built from 3, 1, and 1 primitives respectively.")
    print(" - '*' (first): A set of d-type polarization functions is added to heavy (non-hydrogen) atoms.")
    print(" - '*' (second): A set of p-type polarization functions is added to hydrogen atoms.")
    
    print(f"\nApplying this to {atom_name}, a heavy atom with a 1s core and 2s, 2p valence shells:")

    # Core shell (1s orbital)
    core_s = 6
    print(f"\n1. Core shell (1s): The 1s orbital is described by {core_s} primitives.")

    # Valence s-shell (2s orbital)
    valence_s = 3 + 1 + 1
    print(f"2. Valence s-shell (2s): The 2s orbital is described by 3+1+1 = {valence_s} primitives.")
    
    # Valence p-shells (2p orbitals)
    # There are three 2p orbitals (px, py, pz), and each gets the '311' treatment.
    valence_p_per_orbital = 3 + 1 + 1
    num_p_orbitals = 3
    valence_p_total = num_p_orbitals * valence_p_per_orbital
    print(f"3. Valence p-shells (2p): Each of the three 2p orbitals is described by 3+1+1 = {valence_p_per_orbital} primitives.")
    print(f"   This gives a total of 3 * {valence_p_per_orbital} = {valence_p_total} primitives for the 2p shells.")

    # Polarization functions (d-type for heavy atoms like Carbon)
    # This corresponds to the first '*' in 6-311G**. A standard set has 5 spherical d-functions.
    polarization_d = 5
    print(f"4. Polarization functions ('*'): For a heavy atom, one set of {polarization_d} d-functions is added.")

    # Total calculation
    total_primitives = core_s + valence_s + valence_p_total + polarization_d
    print("\n---------------------------------------------------------------------------------")
    print("The total number of primitive Gaussians for Carbon is the sum of these parts:")
    
    # Final equation as requested
    print(f"Equation: {core_s} (1s) + {valence_s} (2s) + {valence_p_total} (2p) + {polarization_d} (d-pol) = {total_primitives}")
    
    print(f"\nTherefore, a Carbon atom is described by {total_primitives} primitive Gaussians in the 6-311G** basis set.")

    print("\n---------------------------------------------------------------------------------")
    print("For context, for a Hydrogen (H) atom:")
    print(" - No core electrons. The 1s orbital is treated as valence.")
    print(" - Valence 1s orbital ('311' split): 3 + 1 + 1 = 5 primitives.")
    print(" - Polarization (second '*' applies): one set of 3 p-functions = 3 primitives.")
    print(" - Total for H = 5 + 3 = 8 primitives.")

if __name__ == '__main__':
    calculate_basis_set_primitives()
<<<31>>>