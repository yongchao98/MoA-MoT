def solve_photon_entanglement():
    """
    Calculates the polarization of a companion photon based on the conservation of angular momentum.
    """
    # Step 1: Assign numerical values for polarization (helicity).
    right_handed = 1
    left_handed = -1

    # Step 2: The total angular momentum of the system must be conserved.
    # The atom starts in a 1S0 state, which has total angular momentum J=0.
    # Therefore, the sum of the helicities of the two photons must be 0.
    total_angular_momentum = 0

    # Step 3: We measure one photon and find it is right-handed.
    helicity_photon1 = right_handed

    # Step 4: Use the conservation equation to find the helicity of the second photon.
    # The equation is: total_angular_momentum = helicity_photon1 + helicity_photon2
    # Solving for helicity_photon2:
    helicity_photon2 = total_angular_momentum - helicity_photon1
    
    print("Principle: Conservation of Angular Momentum")
    print(f"Initial angular momentum of the atom: {total_angular_momentum}")
    print("Therefore, the total angular momentum of the two emitted photons must also be zero.")
    print("-" * 30)
    print("Representing polarization as helicity:")
    print(f"Right-handed = +{right_handed}")
    print(f"Left-handed = {left_handed}")
    print("-" * 30)
    print("Calculation:")
    print("The conservation equation is: (Helicity of Photon 1) + (Helicity of Photon 2) = Total Angular Momentum")
    print(f"The final equation with known values is: {helicity_photon1} + (Helicity of Photon 2) = {total_angular_momentum}")
    
    print("\nSolving for the companion photon:")
    print(f"Helicity of Photon 2 = {total_angular_momentum} - {helicity_photon1}")
    print(f"Helicity of Photon 2 = {helicity_photon2}")

    # Step 5: Determine the polarization from the resulting helicity.
    if helicity_photon2 == left_handed:
        result = "Left-handed"
    elif helicity_photon2 == right_handed:
        result = "Right-handed"
    else:
        result = "Undetermined"
        
    print(f"\nA helicity of {helicity_photon2} corresponds to a {result} polarization.")

solve_photon_entanglement()