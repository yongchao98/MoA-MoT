def display_density_profile_solution():
    """
    This function prints the formatted analytical solution for the specified
    shock wave density profile.
    """

    # Define the numerical coefficients of the final analytical solution.
    # The equation is ρ/ρ₀ = c1 + c2 / sqrt(c3 + c4 * exp(c5/c6 * x/L))
    c1 = 1
    c2 = 1
    c3 = 1
    c4 = 3
    c5 = -4
    c6 = 3

    # Print the equation in a readable, formatted way
    print("The analytical solution for the density profile is:")
    print()
    print(f"  ρ(x)          {c2}")
    print(f"------- = {c1} + --------------------------------------")
    print("   ρ₀          _________________________________")
    print(f"              ╱        ( {c5}   x )")
    print(f"             ╱  {c3} + {c4} exp( --- * - )")
    print(f"            ╱          ( {c6}   L )")
    print("          ╲╱")
    print()
    print("Where:")
    print("  ρ(x) is the density at position x.")
    print("  ρ₀   is the ambient (pre-shock) density.")
    print("  x    is the position coordinate.")
    print("  L    is the ambient conductive length scale L = κ/(ρ₀⋅M⋅c₀⋅Cᵥ).")
    print("  exp  is the exponential function.")

if __name__ == "__main__":
    display_density_profile_solution()
