def print_fermionic_partition_function():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time representation using Feynman's path integral formalism.
    """

    # --- The Main Formula ---
    print("The formula for the fermionic partition function (Z) is:")
    print("Z = \u222B D\u03C8\u0304 D\u03C8 exp(-S[\u03C8\u0304, \u03C8])")
    print("\n")

    # --- Explanation of Terms ---
    print("Where:")
    print(" - Z: The partition function.")
    print(" - \u03C8(x, \u03C4), \u03C8\u0304(x, \u03C4): Anti-commuting Grassmann fields representing the fermions.")
    print(" - \u222B D\u03C8\u0304 D\u03C8: The functional integral over all possible field configurations.")
    print("\n")

    # --- Boundary Conditions ---
    print("The integral is performed under anti-periodic boundary conditions in imaginary time \u03C4:")
    print("\u03C8(x, \u03C4 = 0) = - \u03C8(x, \u03C4 = \u03B2)")
    print("\n")

    # --- The Euclidean Action S ---
    print("And S[\u03C8\u0304, \u03C8] is the Euclidean action of the system:")
    # Note: Explicitly showing the numbers 0 and β as requested.
    print("S[\u03C8\u0304, \u03C8] = \u222B from 0 to \u03B2 d\u03C4 \u222B d\u1d48x [ \u03C8\u0304(x, \u03C4)(\u2202_\u03C4 + H\u2080)\u03C8(x, \u03C4) + L_int ]")
    print("\n")

    # --- Explanation of Action Components ---
    print("Components of the action S:")
    print(" - \u03C4: Imaginary time variable.")
    print(" - x: Spatial coordinates.")
    print(" - \u03B2 = 1/(k_B T): Inverse temperature, with k_B as Boltzmann's constant and T as temperature.")
    print(" - \u2202_\u03C4: Partial derivative with respect to imaginary time.")
    print(" - H\u2080: The single-particle Hamiltonian operator (e.g., -(\u210F\u00B2/2m)\u2207\u00B2 - \u03BC).")
    print(" - L_int: The interaction part of the Lagrangian density, if any.")

# Execute the function to print the formula
print_fermionic_partition_function()