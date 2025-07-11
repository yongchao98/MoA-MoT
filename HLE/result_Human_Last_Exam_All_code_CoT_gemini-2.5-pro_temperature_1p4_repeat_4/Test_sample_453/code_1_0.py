import sympy

def calculate_energy_shift():
    """
    Calculates the second-order relativistic energy shift for a hydrogen atom
    in the n=3, l=2 state.
    """
    # Quantum numbers for the state
    n = 3
    l = 2
    Z = 1 # For hydrogen atom

    # The calculation simplifies for n = l + 1, as S1 becomes zero.
    # n - l - 1 = 3 - 2 - 1 = 0
    # So the sum S1(n,l) is zero.

    # Calculate the term inside the brackets from the general formula
    term1_denom_factor = l * (l + 1)
    term1 = sympy.Rational(1, 4 * term1_denom_factor**2)

    term2_denom_factor = n * l * (l + 1)
    term2 = sympy.Rational(3, 4 * term2_denom_factor)

    bracket_term = term1 - term2

    # Calculate the overall prefactor
    prefactor = -sympy.Rational(1, n**3)

    # Total numerical coefficient
    coefficient = prefactor * bracket_term

    # We express the final answer in terms of fundamental constants.
    # The full energy shift is Coefficient * m * c**2 * alpha**6
    # alpha = e**2 / (4 * pi * epsilon_0 * hbar * c)

    num = coefficient.p
    den = coefficient.q
    
    print("The second-order energy shift is given by the formula:")
    print("ΔE = C * m * c^2 * α^6")
    print("where C is a numerical coefficient, m is the electron mass,")
    print("c is the speed of light, and α is the fine-structure constant.")
    print("\nCalculation of the coefficient C:")
    print(f"For n={n} and l={l}:")
    print(f"C = - (1/n^3) * [ 1/(4*l(l+1))^2 - 3/(4*n*l(l+1)) ]")
    print(f"C = - (1/{n**3}) * [ 1/(4*{l*(l+1)})^2 - 3/(4*{n}*{l*(l+1)}) ]")
    print(f"C = - (1/{n**3}) * [ 1/{4 * (l*(l+1))**2} - 3/{4*n*l*(l+1)} ]")
    print(f"C = - (1/27) * [ 1/144 - 1/24 ]")
    print(f"C = - (1/27) * [ (1 - 6) / 144 ]")
    print(f"C = - (1/27) * [ -5 / 144 ]")
    print(f"C = {num}/{den}")

    print("\nFinal result expressed in terms of fundamental constants:")
    print(f"ΔE = ({num}/{den}) * m * c^2 * (e^2 / (4 * π * ε_0 * ħ * c))^6")


calculate_energy_shift()