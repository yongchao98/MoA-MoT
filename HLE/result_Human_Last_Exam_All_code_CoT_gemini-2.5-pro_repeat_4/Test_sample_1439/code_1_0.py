import sympy

def find_leading_order_contribution():
    """
    This function determines the order in the coupling constant 'u' at which the
    critical exponent ν receives its first non-vanishing contribution.

    Background:
    1. In the RG analysis of φ⁴ theory, the critical exponent ν is given by:
       1/ν = 2 - γ_r(u*)
       where u* is the coupling at the non-trivial fixed point and γ_r(u) is an RG function.

    2. The function γ_r(u) is calculated as a power series in u via loop expansion.
       γ_r(u) = c1*u + c2*u² + c3*u³ + ...
       The coefficients c1, c2, etc., are derived from 1-loop, 2-loop, etc., diagrams.

    3. The first non-vanishing contribution comes from the 1-loop diagram. This
       diagram has one vertex, which contributes one power of the coupling 'u'.
       Therefore, the coefficient c1 is non-zero, and the expansion starts at order u¹.

    This script formalizes finding this leading order.
    """
    # Define the coupling constant 'u' as a symbolic variable.
    u = sympy.Symbol('u')

    # Represent the non-zero one-loop coefficient as a symbolic constant 'c1'.
    c1 = sympy.Symbol('c1')

    # The perturbative expansion for the RG function γ_r(u) starts at order u^1.
    # We represent this as a series in u, where O(u**2) stands for terms of order u^2 and higher.
    gamma_r_series = c1 * u + sympy.O(u**2)

    # The "initial non-vanishing contribution" corresponds to the leading term of this series.
    leading_term = gamma_r_series.as_leading_term(u)

    # The order of this contribution is the power of 'u' in the leading term.
    # We use sympy's degree function to extract this power.
    order = sympy.degree(leading_term, u)

    print(f"The RG function γ_r(u) is expanded as a series in the coupling constant u.")
    print(f"The series begins: γ_r(u) = {gamma_r_series}")
    print(f"The first non-vanishing contribution is the leading term of this series: {leading_term}")
    print(f"The order of this contribution is the exponent of 'u' in this term.")
    print(f"Order = {order}")

# Execute the function to find and print the result.
find_leading_order_contribution()