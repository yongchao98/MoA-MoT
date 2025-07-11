import math

def calculate_h2_plus_2s_overlap(R, zeta):
    """
    Calculates the overlap integral S for two 2s orbitals in H2+
    for a given internuclear distance R and effective nuclear charge zeta.

    The analytical formula for S(R, ζ) has been derived by integrating the
    product of the two wavefunctions over all space.

    Args:
        R (float): The internuclear distance in atomic units (bohrs).
        zeta (float): The effective nuclear charge.

    Returns:
        float: The value of the overlap integral S.
    """
    # Let a = ζR/2 to simplify the expression
    a = zeta * R / 2.0
    
    # The derived analytical formula is S = e^(-a) * [1 + a + a^2/3 + a^4/15]
    # which expands to the formula in terms of R and ζ shown in the printout below.
    
    polynomial_term = 1.0 + a + (a**2 / 3.0) + (a**4 / 15.0)
    S_value = math.exp(-a) * polynomial_term
    return S_value

def print_overlap_integral_formula():
    """
    Prints the derived analytical formula for the 2s-2s overlap integral,
    showing each numerical coefficient explicitly.
    """
    print("The analytical expression for the overlap integral S between two hydrogenic 2s orbitals is:")
    
    # Define the coefficients of the polynomial part of the formula in terms of (ζ*R).
    # These coefficients result from the multi-step integration in elliptical coordinates.
    c0 = 1.0
    c1 = 1.0 / 2.0
    c2 = 1.0 / 12.0
    c3 = 0.0  # Note: The cubic term's coefficient is zero.
    c4 = 1.0 / 240.0

    print("\nS(R, ζ) = exp(-ζ*R/2) * [ C0 + C1*(ζ*R) + C2*(ζ*R)² + C3*(ζ*R)³ + C4*(ζ*R)⁴ ]")
    
    print("\nThe specific numerical coefficients are:")
    print(f"The coefficient for the (ζ*R)⁰ term is C0 = {c0}")
    print(f"The coefficient for the (ζ*R)¹ term is C1 = {c1}")
    print(f"The coefficient for the (ζ*R)² term is C2 = {c2:.8f} (which is 1/12)")
    print(f"The coefficient for the (ζ*R)³ term is C3 = {c3}")
    print(f"The coefficient for the (ζ*R)⁴ term is C4 = {c4:.8f} (which is 1/240)")

    print("\nTherefore, the final equation with each number explicitly shown is:")
    # We will print the equation using the numerical coefficients
    print(f"S(R, ζ) = exp(-ζ*R/2) * ( {c0} + {c1}*(ζ*R) + {c2:.8f}*(ζ*R)² + {c4:.8f}*(ζ*R)⁴ )")

# --- Main execution ---
# Print the formula as the primary task.
print_overlap_integral_formula()

# As a functional example, let's calculate the overlap for H2+ near its
# equilibrium distance R = 2.0 bohr and with an optimized effective charge ζ ≈ 1.24.
R_example = 2.0  # bohr
zeta_example = 1.24 
S_value = calculate_h2_plus_2s_overlap(R_example, zeta_example)

print(f"\n\n--- Example Calculation ---")
print(f"For an internuclear distance R = {R_example} bohr and an effective charge ζ = {zeta_example}:")
print(f"The calculated overlap integral S is: {S_value:.8f}")
