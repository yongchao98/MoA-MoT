def solve_magnetic_shielding():
    """
    This function prints the derived results for the permeability of a
    cylindrical magnetic shell and the resulting magnetic field inside it.
    """

    # --- Part 1: Required Permeability ---
    print("This problem asks for the properties of a cylindrical shell that perfectly")
    print("shields its exterior from any distortion due to its presence in a uniform magnetic field.")
    print("\nBased on the derivation from magnetostatic boundary conditions, the required")
    print("permeability 'mu' of the shell material is found to be:")

    # Using string variables to represent the physical constants and variables
    mu_0 = "μ₀"
    mu = "μ"

    print("\nFinal Equation for Permeability:")
    print(f"  {mu} = -{mu_0}")
    print("\nThis corresponds to a relative permeability of μ_r = -1. Such a material is")
    print("an idealized 'metamaterial' which is not found in nature for static fields but is")
    print("the required value to achieve perfect external shielding.")

    print("\n" + "="*60 + "\n")

    # --- Part 2: Magnetic Field in the Interior ---
    print("With this permeability, the magnetic field 'H_int' in the interior region (ρ < R₁) is:")

    # Defining variables for the equation
    H0 = "H₀"
    R1 = "R₁"
    R2 = "R₂"
    x_hat = "x̂"  # Unicode for x-hat vector symbol

    print("\nFinal Equation for Interior Magnetic Field:")
    # Printing the equation components as requested
    print(f"  H_int = ( {H0} ) * ( {R2}² / {R1}² ) * ( {x_hat} )")
    print("\nThis shows that the field inside the shell is uniform and points in the same direction")
    print(f"as the applied external field, but its magnitude is amplified by the factor ({R2}/{R1})².")

solve_magnetic_shielding()