import sympy
from sympy import exp, oo, integrate, pi, pprint, Symbol, Rational

def solve_overlap_integral():
    """
    This function symbolically derives the overlap integral for two 2s orbitals
    in a diatomic molecule using elliptical coordinates and prints the
    step-by-step derivation and the final analytical expression.
    """
    print("Deriving the overlap integral for two 2s orbitals S = \u222B \u03C8_2s(A) * \u03C8_2s(B) d\u03C4\n")
    
    print("Step 1: Define symbols and coordinate system.")
    # Use 'l' for lambda as 'lambda' is a reserved keyword in Python.
    zeta, R = sympy.symbols('zeta R', positive=True, real=True)
    l, m, phi = sympy.symbols('lambda mu phi', real=True)
    
    # For simplification, we introduce rho = zeta * R
    rho = Symbol('rho', positive=True, real=True)

    print("The transformation to elliptical coordinates is:")
    print("r_A = R/2 * (\u03BB + \u03BC)")
    print("r_B = R/2 * (\u03BB - \u03BC)")
    print("d\u03C4 = (R\u00B3/8) * (\u03BB\u00B2 - \u03BC\u00B2) d\u03BB d\u03BC d\u03C6")
    print("-" * 40)

    print("Step 2: Construct the integrand.")
    # The normalized 2s wave function (ψ_2s) is proportional to (2 - ζr) * exp(-ζr / 2).
    # The product ψ_2s(A)ψ_2s(B) in elliptical coordinates, in terms of ρ = ζR, becomes:
    # [4 - 2ρλ + (ρ²/4)(λ²-μ²)] * exp(-ρλ/2)
    poly_term = 4 - 2*rho*l + (rho**2 / 4)*(l**2 - m**2)
    exp_term = exp(-rho*l / 2)
    volume_factor_part = (l**2 - m**2)
    
    integrand = poly_term * exp_term * volume_factor_part
    print("The core part of the integrand (a function of \u03BB, \u03BC, \u03C1) is:")
    pprint(integrand)
    print("-" * 40)

    print("Step 3: Integrate over \u03C6 from 0 to 2\u03C0.")
    integral_phi = integrate(1, (phi, 0, 2*pi))
    print(f"\u222B d\u03C6 from 0 to 2\u03C0 = {integral_phi}")
    print("-" * 40)

    print("Step 4: Integrate over \u03BC from -1 to 1.")
    integral_mu = integrate(integrand, (m, -1, 1))
    integral_mu = sympy.simplify(integral_mu)
    print("Result after integrating with respect to \u03BC:")
    pprint(integral_mu)
    print("-" * 40)

    print("Step 5: Integrate over \u03BB from 1 to \u221E.")
    integral_lambda = integrate(integral_mu, (l, 1, oo))
    integral_lambda = sympy.simplify(integral_lambda)
    print("Result after integrating with respect to \u03BB:")
    pprint(integral_lambda)
    print("-" * 40)

    print("Step 6: Combine all factors to get the final overlap integral S(\u03C1).")
    # Constant factor C = (Normalization const)² * (Volume element const) * (phi integral)
    # C = (ζ³/32π) * (R³/8) * (2π) = ζ³R³ / 128 = ρ³ / 128
    constant_factor = rho**3 / 128
    
    S_rho_expanded = sympy.expand(sympy.simplify(constant_factor * integral_lambda))
    
    print("The final expression for S in terms of \u03C1 (where \u03C1 = \u03B6R) is:")
    pprint(S_rho_expanded)
    print("-" * 40)
    
    print("\n--- FINAL ANALYTICAL EXPRESSION ---")
    print("The overlap integral S is a function of \u03C1 = \u03B6R.")
    
    # Deconstruct the expression for formatted printing
    c0 = S_rho_expanded.coeff(rho, 0)
    c1 = S_rho_expanded.coeff(rho, 1)
    c2 = S_rho_expanded.coeff(rho, 2)
    c4 = S_rho_expanded.coeff(rho, 4)
    
    print("The equation is of the form: S(\u03C1) = exp(-\u03C1/2) * [c0 + c1*\u03C1 + c2*\u03C1\u00B2 + c4*\u03C1\u2074]")

    print("\nThus, the final equation with each number explicitly shown is:")
    print(f"S(\u03C1) = exp(-\u03C1 / {2}) * ( {c0} + ({c1.p}/{c1.q})*\u03C1 + ({c2.p}/{c2.q})*\u03C1\u00B2 + ({c4.p}/{c4.q})*\u03C1\u2074 )")
    print("where \u03C1 = \u03B6 * R")

if __name__ == '__main__':
    solve_overlap_integral()
<<<S(rho) = exp(-rho / 2) * ( 1 + (1/2)*rho + (1/12)*rho**2 + (1/240)*rho**4 )>>>