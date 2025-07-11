def analyze_ammonia_tunneling():
    """
    Analyzes whether an ammonia molecule with exotic spin-0 hydrogen nuclei
    would exhibit quantum tunneling.
    """

    print("--- Step 1: Defining the Particles ---")
    print("Ordinary hydrogen nuclei (protons) are fermions.")
    ordinary_spin = 1/2
    print(f"They have nuclear spin = {ordinary_spin}.")
    print("The exotic hydrogen nuclei are defined as bosons.")
    exotic_spin = 0
    print(f"They have nuclear spin = {exotic_spin}.")
    print("\n")

    print("--- Step 2: The Pauli Exclusion Principle ---")
    print("The Pauli principle dictates the symmetry of the total wavefunction (Ψ_total) when identical particles are exchanged.")
    print(" - For a system of identical FERMIONS (like ordinary protons), Ψ_total must be ANTISYMMETRIC.")
    print(" - For a system of identical BOSONS (like exotic H-nuclei), Ψ_total must be SYMMETRIC.")
    print("\n")

    print("--- Step 3: Analyzing the Wavefunction for Exotic Ammonia (N(H_exotic)3) ---")
    print("The total wavefunction can be approximated as a product:")
    print("Ψ_total = Ψ_nuclear_spin * Ψ_rotational * Ψ_vibrational")
    print("The Pauli principle for the three identical exotic H-nuclei (bosons) requires Ψ_total to be SYMMETRIC.\n")
    
    print("--- Step 4: Analyzing the Symmetry of Each Part ---")
    
    # Nuclear Spin Part
    print("1. Ψ_nuclear_spin:")
    print("Since the exotic nuclei have spin 0, the total nuclear spin is 0.")
    print("There is only one possible nuclear spin state, and it is always SYMMETRIC under particle exchange.")
    print("\n")

    # The consequence of the spin part being symmetric
    print("Because Ψ_total must be SYMMETRIC, and Ψ_nuclear_spin is already SYMMETRIC, the product of the remaining parts must also be SYMMETRIC:")
    print("Symmetry(Ψ_rotational * Ψ_vibrational) must be SYMMETRIC.")
    print("\n")

    # Vibrational Inversion Part
    print("2. Ψ_vibrational (for the inversion motion):")
    print("Tunneling of the Nitrogen atom through the plane of the hydrogens splits the ground state into a doublet:")
    print(" - A lower energy state with a SYMMETRIC vibrational wavefunction (let's call it Vib_S).")
    print(" - A higher energy state with an ANTISYMMETRIC vibrational wavefunction (let's call it Vib_A).")
    print("The existence of both of these states IS the phenomenon of tunneling.\n")

    # Rotational Part
    print("3. Ψ_rotational:")
    print("The rotational wavefunctions of the molecule can also be classified by their symmetry under exchange of the three H-nuclei.")
    print("The molecule can exist in rotational states that are SYMMETRIC (Rot_S) or ANTISYMMETRIC (Rot_A) with respect to this exchange.")
    print("\n")

    print("--- Step 5: Finding Pauli-Allowed Combinations ---")
    print("The product Ψ_rotational * Ψ_vibrational must be SYMMETRIC.")
    print("We can achieve this in two ways:")
    print(" - Combination 1: (Symmetric Rotational State) * (Symmetric Vibrational State)")
    print("   Rot_S * Vib_S  -->  SYMMETRIC. This is ALLOWED.")
    print("\n")
    print(" - Combination 2: (Antisymmetric Rotational State) * (Antisymmetric Vibrational State)")
    print("   Rot_A * Vib_A  -->  SYMMETRIC. This is also ALLOWED.")
    print("\n")

    print("--- Step 6: Conclusion ---")
    print("Because allowed rotational states exist for BOTH the symmetric and antisymmetric vibrational states of the tunneling doublet, BOTH energy levels are physically realizable.")
    print("The molecule can populate the lower (symmetric) tunneling state when it is in a symmetric rotational state, and it can populate the upper (antisymmetric) tunneling state when it is in an antisymmetric rotational state.")
    print("Since both states of the tunneling doublet can exist, the energy splitting exists.")
    print("\n")
    print("Therefore, the ammonia molecule with exotic hydrogens would still exhibit tunneling.")

if __name__ == '__main__':
    analyze_ammonia_tunneling()
