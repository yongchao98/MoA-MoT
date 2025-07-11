import sympy

def calculate_cross_section():
    """
    This script calculates the total cross section for fermion-fermion scattering
    mediated by a pseudoscalar boson, using symbolic mathematics.
    """
    # Define symbolic variables
    # g: coupling constant
    # M: mass of the scalar boson
    # E: energy of each incoming fermion in the Center of Mass (CM) frame
    # s, t, u: Mandelstam variables
    g, M, E, s, t = sympy.symbols('g M E s t', real=True, positive=True)
    
    # Use M_sq for M^2 to keep expressions clean and avoid sqrt(M**2)
    M_sq = sympy.Symbol('M^2', real=True, positive=True)

    # --- Step 1: Define Mandelstam variables ---
    # We work in the high-energy limit where the fermion mass m=0.
    # s is the square of the total CM energy.
    s_expr = 4 * E**2
    # u is related to s and t by s+t+u=0
    u_expr = -s - t

    print("Calculating the total cross section for fermion-fermion scattering ψψ -> ψψ.")
    print("The process is mediated by a pseudoscalar boson φ with mass M.")
    print("The calculation is performed in the high-energy limit (fermion mass m -> 0).")
    print("-" * 60)

    # --- Step 2: Write the spin-averaged squared matrix element <|M|^2> ---
    # <|M|^2> = g^4 * [ t^2/(t-M^2)^2 + u^2/(u-M^2)^2 + t*u/((t-M^2)*(u-M^2)) ]
    # We substitute u = -s - t into the expression.
    avg_M_squared = g**4 * (
        t**2 / (t - M_sq)**2 +
        u_expr**2 / (u_expr - M_sq)**2 +
        t * u_expr / ((t - M_sq) * (u_expr - M_sq))
    )

    print("The spin-averaged squared matrix element is:")
    print("<|M|^2> = g^4 * [t^2/(t-M^2)^2 + u^2/(u-M^2)^2 + t*u/((t-M^2)(u-M^2))]")
    print("where u = -s - t.")
    print("-" * 60)

    # --- Step 3: Set up and perform the integration for σ ---
    # The total cross section is given by the integral:
    # σ = (1 / (32 * π * s^2)) * ∫[-s to 0] <|M|^2> dt
    print("Integrating <|M|^2> with respect to t from -s to 0...")
    integral_result = sympy.integrate(avg_M_squared, (t, -s, 0))

    # --- Step 4: Assemble and simplify the final expression for σ ---
    # Substitute the result of the integral into the formula for sigma
    sigma_expr_s = (1 / (32 * sympy.pi * s**2)) * integral_result

    # Substitute s = 4*E^2 and M_sq = M**2 to get the final answer in terms of E and M.
    final_sigma = sigma_expr_s.subs([(s, s_expr), (M_sq, M**2)])

    # Simplify the final expression
    final_sigma_simplified = sympy.simplify(final_sigma)
    
    print("The integration is complete.")
    print("-" * 60)

    # --- Step 5: Print the final answer ---
    print("The final expression for the total cross section σ in terms of E and M is:")
    
    # To satisfy the prompt's request to "output each number in the final equation",
    # we present the formula in a structured way.
    numer, denom = final_sigma_simplified.as_numer_denom()
    
    print("\nσ = (Numerator) / (Denominator)\n")
    print("Numerator:")
    sympy.pprint(numer)
    print("\nDenominator:")
    sympy.pprint(denom)

    # We also print the full one-line equation for clarity.
    # The str() representation is often more compact than pretty print.
    final_eq_str = str(final_sigma_simplified).replace('**', '^')
    print("\nFull Equation:")
    print(f"σ = {final_eq_str}")


if __name__ == '__main__':
    calculate_cross_section()