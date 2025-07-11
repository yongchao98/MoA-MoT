import numpy as np

def calculate_2s_2s_overlap():
    """
    Calculates the overlap integral S for two hydrogenic 2s orbitals and prints the
    derivation formula and a numerical example.
    """
    # Define example parameters for the calculation.
    # The problem asks for an expression in terms of R and zeta. The code will
    # first display this general formula and then compute a specific case.
    # We will use zeta=1.0 for a pure hydrogenic orbital and an example distance R.
    zeta = 1.0
    # A typical internuclear distance for diatomic molecules is around 2.0 a.u. (Bohr radii)
    R = 2.0

    print("This program calculates the overlap integral S for two hydrogenic 2s orbitals.")
    print("The analytical expression is derived using elliptical coordinates.")
    print("-" * 70)

    # Print the analytical formula with its coefficients
    print("The analytical expression for the overlap integral S as a function of ρ = ζR is:")
    # The prompt requests to output each number in the final equation.
    # The polynomial is P(ρ) = c0 + c1*ρ + c2*ρ^2 + c3*ρ^3 + c4*ρ^4
    # The coefficients are c0=1, c1=1/2, c2=1/12, c3=0, c4=1/240
    print("S(ρ) = exp(-ρ/2) * (1 + (1/2)ρ + (1/12)ρ² + (1/240)ρ⁴)")
    print("\nThe numbers (coefficients for ρ⁰, ρ¹, ρ², ρ³, ρ⁴) in the polynomial part are:")
    print("c0 = 1")
    print("c1 = 1/2  = 0.5")
    print("c2 = 1/12 ≈ 0.0833")
    print("c3 = 0")
    print("c4 = 1/240 ≈ 0.00417")

    print("-" * 70)
    print(f"Calculating for the example case: ζ = {zeta}, R = {R} a.u.")

    # Calculate rho
    rho = zeta * R
    print(f"\nFirst, calculate ρ = ζ * R = {zeta:.2f} * {R:.2f} = {rho:.4f}")

    # Substitute rho into the formula and show the calculation step-by-step
    term_exp = np.exp(-rho / 2.0)
    term0 = 1.0
    term1 = rho / 2.0
    term2 = rho**2 / 12.0
    term4 = rho**4 / 240.0
    
    poly_sum = term0 + term1 + term2 + term4
    S = term_exp * poly_sum

    print("\nNext, substitute the value of ρ into the equation:")
    print(f"S = exp(-{rho:.2f}/2) * (1 + {rho:.2f}/2 + {rho:.2f}²/12 + {rho:.2f}⁴/240)")
    print(f"S = exp({-rho/2.0:.4f}) * (1 + {term1:.4f} + {rho**2:.4f}/12 + {rho**4:.4f}/240)")
    print(f"S = {term_exp:.6f} * (1 + {term1:.4f} + {term2:.4f} + {term4:.4f})")
    print(f"S = {term_exp:.6f} * ({poly_sum:.6f})")
    
    # Print the final numerical result
    print("-" * 70)
    print(f"The final calculated overlap integral is: S = {S:.6f}")


# Execute the main function
calculate_2s_2s_overlap()