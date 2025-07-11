def solve_photon_entanglement():
    """
    Calculates the polarization of a companion photon based on the conservation of angular momentum.
    """
    # Define the mapping from spin projection to polarization name
    polarization_map = {
        1: "Right-handed",
        -1: "Left-handed"
    }

    # Initial and final total angular momentum (J) of the atom are both 0 (from the 1S0 state)
    J_initial_atom = 0
    J_final_atom = 0

    # The first photon is measured to be right-handed, which corresponds to a spin projection of +1
    spin_photon1 = 1
    polarization_photon1 = polarization_map[spin_photon1]

    # The law of conservation of angular momentum states:
    # J_initial_atom = J_final_atom + spin_photon1 + spin_photon2
    # So, spin_photon2 = J_initial_atom - J_final_atom - spin_photon1
    spin_photon2 = J_initial_atom - J_final_atom - spin_photon1

    # Determine the polarization of the second photon from its spin
    polarization_photon2 = polarization_map[spin_photon2]

    print("This problem is governed by the conservation of angular momentum.")
    print(f"The initial state of the atom has total angular momentum J = {J_initial_atom}.")
    print(f"The final state of the atom has total angular momentum J = {J_final_atom}.")
    print(f"The first photon is measured as {polarization_photon1}, which corresponds to a spin of {spin_photon1}.")
    print("\nTo conserve total angular momentum, the following equation must be true:")
    print("Initial_J = Final_J + Spin_Photon1 + Spin_Photon2")
    # Outputting each number in the final equation
    print(f"   {J_initial_atom}    =    {J_final_atom}    +    ({spin_photon1})    +    ({spin_photon2})")
    
    print(f"\nTherefore, the second photon must have a spin of {spin_photon2}.")
    print(f"A spin of {spin_photon2} corresponds to a polarization that is: {polarization_photon2}.")

solve_photon_entanglement()