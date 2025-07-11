import sympy
from sympy import exp, integrate, oo, pi, sqrt, pretty_print

def solve_overlap_integral():
    """
    Calculates the overlap integral for two 2s orbitals in H2+
    using symbolic mathematics.
    """
    # Step 1 & 2: Define symbols and set up the integral components
    print("1. Defining symbols and integral components...")
    R, zeta = sympy.symbols('R zeta', positive=True, real=True)
    lam, mu, phi = sympy.symbols('lambda mu phi', real=True)

    # The dimensionless variable rho = zeta * R
    rho = zeta * R

    # Normalization constant for a 2s orbital squared
    # N = (1/(4*sqrt(2*pi))) * zeta**(3/2) -> N**2 = zeta**3 / (32*pi)
    norm_const_sq = zeta**3 / (32 * pi)

    # The product of the spatial parts of the two wavefunctions in elliptical coordinates
    # (2 - zeta*r_A) * (2 - zeta*r_B) = 4 - 2*zeta*(r_A+r_B) + zeta**2*r_A*r_B
    # In elliptical coords: r_A+r_B = R*lam and r_A*r_B = (R**2/4)*(lam**2 - mu**2)
    # This becomes: 4 - 2*rho*lam + (rho**2/4)*(lam**2 - mu**2)
    wavefunc_part = 4 - 2 * rho * lam + (rho**2 / 4) * (lam**2 - mu**2)

    # The exponential part of the wavefunction product
    # exp(-zeta*r_A/2) * exp(-zeta*r_B/2) = exp(-zeta*(r_A+r_B)/2) = exp(-rho*lam/2)
    exp_part = exp(-rho * lam / 2)

    # The volume element in elliptical coordinates
    # d_tau = (R**3/8) * (lam**2 - mu**2) d_lam d_mu d_phi
    volume_element = (R**3 / 8) * (lam**2 - mu**2)

    # Combine all parts to form the full integrand
    integrand = norm_const_sq * wavefunc_part * exp_part * volume_element

    # Step 3: Perform the integration sequentially
    print("\n2. Integrating over phi, mu, and lambda...")

    # Integrate over phi from 0 to 2*pi (azimuthal angle)
    # The integrand is independent of phi, so this is a multiplication by 2*pi
    integral_phi = integrate(integrand, (phi, 0, 2 * pi))

    # Integrate over mu from -1 to 1
    integral_mu = integrate(integral_phi, (mu, -1, 1))

    # Integrate over lambda from 1 to infinity
    S = integrate(integral_mu, (lam, 1, oo))

    # Step 4: Simplify and display the final result
    print("\n3. Final analytical expression for the overlap integral S(rho):")
    # The result is in terms of R and zeta, which can be simplified to be in terms of rho
    S_simplified = sympy.simplify(S)
    
    # To make the expression cleaner, let's substitute R = rho/zeta
    S_final = S_simplified.subs(R, rho/zeta).simplify()

    # Print the final equation
    print("\nThe overlap integral S as a function of rho (where rho = zeta * R) is:")
    pretty_print(sympy.Eq(sympy.Symbol('S'), S_final))

    # Step 5: Output the coefficients of the polynomial part of the expression
    print("\nThe expression is of the form: exp(-rho/2) * P(rho)")
    print("The coefficients for the polynomial P(rho) are:")
    
    # The expression is exp(-rho/2) * (polynomial in rho)
    # We extract the polynomial part
    poly_part = sympy.expand(S_final / exp(-rho/2))
    
    # Get coefficients for each power of rho
    coeff_rho4 = poly_part.coeff(rho, 4)
    coeff_rho3 = poly_part.coeff(rho, 3)
    coeff_rho2 = poly_part.coeff(rho, 2)
    coeff_rho1 = poly_part.coeff(rho, 1)
    coeff_rho0 = poly_part.coeff(rho, 0)

    print(f"  - Coefficient of rho^4: {coeff_rho4} = 1/240")
    print(f"  - Coefficient of rho^3: {coeff_rho3}")
    print(f"  - Coefficient of rho^2: {coeff_rho2} = 1/12")
    print(f"  - Coefficient of rho^1: {coeff_rho1} = 1/2")
    print(f"  - Constant term (rho^0): {coeff_rho0}")

solve_overlap_integral()