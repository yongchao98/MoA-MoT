def solve_photon_entanglement():
    """
    This function calculates the polarization of a companion photon
    based on the principle of conservation of angular momentum.
    """
    # In the S0 state, the total angular momentum (J) of the atom is 0.
    J_atom_initial = 0
    J_atom_final = 0

    # By conservation of angular momentum, the total angular momentum of the
    # emitted photons must be the difference, which is 0.
    J_total_photons = J_atom_initial - J_atom_final
    
    # Let's represent the spin of photons:
    # Right-handed polarization -> spin = +1
    # Left-handed polarization -> spin = -1
    polarization_map = {
        1: "Right-handed",
        -1: "Left-handed"
    }

    # The first photon is measured to be right-handed.
    spin_photon1 = 1
    measured_polarization_photon1 = polarization_map[spin_photon1]
    
    # The total spin of the two photons must sum to J_total_photons (which is 0).
    # J_total_photons = spin_photon1 + spin_photon2
    # Therefore, spin_photon2 = J_total_photons - spin_photon1
    spin_photon2 = J_total_photons - spin_photon1
    
    # Determine the polarization of the second photon.
    companion_polarization = polarization_map[spin_photon2]

    print("Step 1: Define initial and final atomic angular momentum.")
    print(f"Initial atomic angular momentum (J_initial): {J_atom_initial}")
    print(f"Final atomic angular momentum (J_final): {J_atom_final}\n")

    print("Step 2: Apply the conservation of angular momentum.")
    print("The total angular momentum of the two photons (J_photons) must be zero.")
    print(f"J_photons = J_initial - J_final = {J_total_photons}\n")

    print("Step 3: A measurement is made on the first photon.")
    print(f"The first photon is measured as {measured_polarization_photon1}, which corresponds to a spin of {spin_photon1}.\n")

    print("Step 4: Calculate the spin of the companion photon.")
    print("The conservation equation is: J_photons = spin_photon1 + spin_photon2")
    # The final equation output as requested
    print(f"The final equation is: {J_total_photons} = {spin_photon1} + ({spin_photon2})")
    
    print(f"\nTherefore, the spin of the second photon must be {spin_photon2}.")
    print(f"This spin value corresponds to a polarization of: {companion_polarization}.\n")

solve_photon_entanglement()