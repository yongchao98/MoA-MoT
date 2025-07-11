def solve_photon_entanglement():
    """
    Solves for the polarization of a companion photon in an entangled pair
    based on the principle of conservation of angular momentum.
    """
    # The atom's initial and final states are J=0.
    # By conservation of angular momentum, the total angular momentum
    # of the two emitted photons must be zero.
    total_angular_momentum = 0

    # We represent circular polarization with helicity values:
    # Right-handed polarization -> Helicity = +1
    # Left-handed polarization  -> Helicity = -1
    polarization_map = {
        1: "Right-handed",
        -1: "Left-handed"
    }

    # One photon is measured to be right-handed.
    helicity_photon_1 = 1
    polarization_photon_1 = polarization_map[helicity_photon_1]

    # The conservation law states:
    # Total Angular Momentum = (Helicity of Photon 1) + (Helicity of Photon 2)
    # 0 = 1 + (Helicity of Photon 2)
    # Therefore, Helicity of Photon 2 = -1
    helicity_photon_2 = total_angular_momentum - helicity_photon_1

    # Find the corresponding polarization of the companion photon.
    polarization_photon_2 = polarization_map[helicity_photon_2]

    # Print the explanation and the final equation
    print("Principle: Conservation of Angular Momentum")
    print(f"Initial Total Angular Momentum (J) = {total_angular_momentum}")
    print(f"Helicity of Measured Photon 1 ({polarization_photon_1}) = {helicity_photon_1}")
    print("\nConservation Equation:")
    print(f"Total Angular Momentum = Helicity of Photon 1 + Helicity of Photon 2")
    # Output the numbers in the final equation as requested
    print(f"{total_angular_momentum} = {helicity_photon_1} + ({helicity_photon_2})")
    print("\nResult:")
    print(f"The calculated helicity of the companion photon is {helicity_photon_2}.")
    print(f"This corresponds to a polarization of: {polarization_photon_2}")

solve_photon_entanglement()
<<<A>>>