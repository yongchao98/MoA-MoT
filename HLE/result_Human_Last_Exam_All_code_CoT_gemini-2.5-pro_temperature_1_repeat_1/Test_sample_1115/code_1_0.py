def solve_photon_polarization():
    """
    Calculates the polarization of a companion photon in an entangled pair
    based on the conservation of angular momentum.
    """
    # Define spin values for polarizations
    # Right-handed polarization corresponds to spin +1
    # Left-handed polarization corresponds to spin -1
    polarization_to_spin = {
        "Right-handed": 1,
        "Left-handed": -1
    }
    spin_to_polarization = {v: k for k, v in polarization_to_spin.items()}

    # Initial and final atomic states have total angular momentum J=0.
    # By conservation of angular momentum, the total spin of the two photons must be 0.
    total_spin = 0

    # One photon is measured to be right-handed.
    measured_photon_polarization = "Right-handed"
    measured_photon_spin = polarization_to_spin[measured_photon_polarization]

    # Calculate the spin of the companion photon using the conservation law.
    # Total Spin = Spin of Photon 1 + Spin of Photon 2
    companion_photon_spin = total_spin - measured_photon_spin

    # Determine the polarization of the companion photon from its spin.
    companion_photon_polarization = spin_to_polarization[companion_photon_spin]

    print("Principle: Conservation of Angular Momentum")
    print("The total spin of the two-photon system must be zero.")
    print("\nThe equation is: Total Spin = (Spin of Photon 1) + (Spin of Photon 2)")
    print(f"Substituting the known values: {total_spin} = ({measured_photon_spin}) + (Spin of Photon 2)")
    print(f"Solving for the companion photon's spin gives: {companion_photon_spin}")
    print(f"\nA spin of {companion_photon_spin} corresponds to {companion_photon_polarization} polarization.")
    print("\nTherefore, the value of the polarization of the companion photon is Left-handed.")


solve_photon_polarization()
<<<A>>>