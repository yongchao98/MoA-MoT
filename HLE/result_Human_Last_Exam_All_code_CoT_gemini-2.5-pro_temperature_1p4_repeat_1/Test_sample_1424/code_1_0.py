import sympy

def display_green_function_dependence():
    """
    This function explains and prints the formula for the bare Green's function,
    highlighting its dependence on single-particle energy eigenvalues.
    """

    # Using sympy to represent the mathematical symbols for a clearer output
    G_0, k, omega, epsilon_k, eta, i = sympy.symbols("G_0 k ω ϵ_k η i")

    # The equation for the bare Green's function
    # G_0(k, ω) = 1 / (ω - ϵ_k + iη)
    # The left-hand side of the equation
    lhs = G_0(k, omega)
    
    # The right-hand side of the equation
    rhs = 1 / (omega - epsilon_k + i*eta)

    print("In the Feynman path integral formalism, the bare Green's function G_0 has a clear and direct dependence on the single-particle energy eigenvalues ϵ_k.")
    print("The relationship is most clearly expressed in frequency (ω) and momentum (k) space.")
    print("\nThe fundamental equation is:")
    
    # We will manually format the print output for clarity as sympy's pretty print can be complex.
    # The equation is: G_0(k, ω) = 1 / (ω - ϵ_k + iη)
    
    equation_lhs = "G_0(k, ω)"
    numerator = "1"
    denominator_part1 = "ω"
    denominator_part2 = "ϵ_k"
    denominator_part3 = "iη"

    print(f"\n  {equation_lhs} = {numerator} / ({denominator_part1} - {denominator_part2} + {denominator_part3})\n")
    
    print("Let's break down each component of this equation:")
    print("-" * 50)
    print(f"Component: {equation_lhs}")
    print("  - Meaning: The bare Green's function for a particle in a state 'k' with frequency 'ω'. It represents the amplitude for the particle to propagate.")
    print("\nComponent: The numerator '1'")
    print(f"  - Meaning: The numerator of the expression is the number 1.")
    print(f"\nComponent: The denominator term '{denominator_part1}'")
    print("  - Meaning: This is 'ω', the frequency variable.")
    print(f"\nComponent: The denominator term '{denominator_part2}'")
    print("  - Meaning: This is 'ϵ_k', the single-particle energy eigenvalue for the state 'k'. This is the core of the dependence.")
    print(f"\nComponent: The denominator term '{denominator_part3}'")
    print("  - Meaning: 'iη' is an infinitesimal imaginary number (where η -> 0+). It is a mathematical device essential for the time-ordering in the Feynman formalism, ensuring causality.")
    print("-" * 50)

    print("\nFunctional Dependence Explained:")
    print(f"The functional dependence of G_0 on ϵ_k is an inverse relationship. Specifically, the Green's function has a pole (the denominator becomes zero) when ω ≈ ϵ_k.")
    print("This pole at the particle's energy is the signature of a stable, non-interacting particle that propagates with its specific energy eigenvalue.")

if __name__ == '__main__':
    display_green_function_dependence()