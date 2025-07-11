import sys

def solve_polarization():
    """
    Calculates the polarization of a companion photon in an entangled pair
    based on the principle of conservation of angular momentum.
    """

    # The initial and final states of the atom have total angular momentum J=0.
    # Therefore, the total angular momentum of the two-photon system must be zero.
    # We consider the projection of spin angular momentum onto the emission axis (z-axis).
    S_z_total = 0

    # Photon 1 is measured to be right-handed. It travels along the +z direction.
    # By definition, its spin is parallel to its momentum.
    # Its spin projection on the z-axis is +1 (in units of ħ).
    S_z_photon1 = 1

    # The conservation of angular momentum gives the equation:
    # S_z_photon1 + S_z_photon2 = S_z_total
    
    # We solve for the spin projection of Photon 2.
    S_z_photon2 = S_z_total - S_z_photon1

    print("--- The Physics ---")
    print("The system is governed by the Conservation of Angular Momentum.")
    print("The initial Calcium atom state has angular momentum J=0.")
    print("Therefore, the total angular momentum of the two emitted photons must also be 0.")
    print("\n--- The Calculation ---")
    print("Let's define the conservation equation for the spin projection (S_z) on the emission axis:")
    print(f"S_z_photon1 + S_z_photon2 = S_z_total")
    print(f"We are given S_z_photon1 = {S_z_photon1} (for a right-handed photon) and S_z_total = {S_z_total}.")
    
    print("\nSolving the equation for the companion photon:")
    # This loop prints out the numbers in the final equation as requested
    equation_str_vars = [str(S_z_photon1), "+ S_z_photon2", "=", str(S_z_total)]
    print(" ".join(equation_str_vars))
    print(f"S_z_photon2 = {S_z_total} - {S_z_photon1}")
    print(f"S_z_photon2 = {S_z_photon2}")

    print("\n--- The Interpretation ---")
    print(f"The companion photon (Photon 2) has a spin projection of {S_z_photon2} on the z-axis.")
    print("However, it travels in the opposite direction (the -z direction).")
    print("Its spin vector points along the -z axis, and its momentum vector also points along the -z axis.")
    print("Since its spin is PARALLEL to its direction of motion, it is, by definition, Right-handed.")

    print("\n--- Final Answer ---")
    print("The value of the polarization of the companion photon is: Right-handed")

solve_polarization()
<<<C>>>