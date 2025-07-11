def solve_photon_polarization():
    """
    Calculates the polarization of a companion photon in an entangled pair
    based on the conservation of angular momentum.
    """
    # The atom's initial and final states (1S0) have total angular momentum J=0.
    J_total = 0

    # We represent the spin projection (polarization) of the photons numerically:
    # Right-handed polarization = +1
    # Left-handed polarization = -1
    # One photon is measured to be right-handed.
    J_photon1 = 1

    # The conservation of angular momentum requires:
    # J_total = J_photon1 + J_photon2
    # We need to find the value of J_photon2.
    J_photon2 = J_total - J_photon1

    # Print the explanation and the final equation with its numbers.
    print("This scenario is governed by the conservation of angular momentum.")
    print(f"The total angular momentum of the system must be conserved, so the sum of the photon spins must be {J_total}.")
    print("The measurement of the first photon found it to be right-handed (spin = +1).")
    print("\nThe conservation equation is: J_total = J_photon1 + J_photon2")
    print("Plugging in the known values, we get the final equation:")
    print(f"{J_total} = {J_photon1} + J_photon2")

    print("\nSolving for the second photon's spin (J_photon2):")
    print(f"J_photon2 = {J_total} - {J_photon1} = {J_photon2}")

    # Interpret the result
    if J_photon2 == -1:
        result_text = "Left-handed"
    elif J_photon2 == 1:
        result_text = "Right-handed"
    else:
        result_text = "Undetermined"
    
    print(f"\nA spin value of {J_photon2} corresponds to {result_text} polarization.")
    print("\nConclusion: The polarization of the companion photon must be Left-handed.")

solve_photon_polarization()