def solve_photon_entanglement():
    """
    Solves the photon entanglement problem based on conservation laws in a J=0 -> J=1 -> J=0 cascade.
    """
    
    # 1. State the physical principle
    print("Principle: In a J=0 -> J=1 -> J=0 atomic cascade, conservation of angular momentum and parity")
    print("requires the two emitted photons to have the same helicity (polarization).")
    print("This means if one is Right-handed, the other must be Right-handed.")
    print("If one is Left-handed, the other must be Left-handed.")
    print("-" * 50)
    
    # 2. Define the measurement performed on the first photon
    measurement_photon_1 = "Right-handed"
    print(f"A measurement is performed on the first photon.")
    print(f"Result: The polarization of photon 1 is '{measurement_photon_1}'.")
    print("-" * 50)

    # 3. Deduce the state of the second photon based on the principle
    if measurement_photon_1 == "Right-handed":
        polarization_photon_2 = "Right-handed"
    elif measurement_photon_1 == "Left-handed":
        polarization_photon_2 = "Left-handed"
    else:
        polarization_photon_2 = "Undetermined by this principle"
        
    print("Deduction: Due to the conservation principle, the polarization of the second photon is determined.")
    
    # 4. Print the final 'equation' showing the result
    print("\nFinal Equation based on Conservation Law:")
    print(f"Polarization(Photon 2) = Polarization(Photon 1)")
    print(f"Polarization(Photon 2) = {measurement_photon_1}")
    print("\nTherefore, the value of the polarization of the companion photon is:")
    print(polarization_photon_2)

solve_photon_entanglement()