import sympy

def mori_tanaka_expression():
    """
    This function symbolically derives and prints the expression for the effective
    elastic moduli (C) in the Mori-Tanaka model for composites.
    """
    # Define the fourth-order tensors as non-commutative symbols, as tensor
    # multiplication is generally not commutative.
    C, Cm, Cf, I, A = sympy.symbols('C Cm Cf I A', commutative=False)

    # Define the volume fractions as standard (commutative) symbols.
    Vf, Vm = sympy.symbols('Vf Vm')

    # Construct the expression for C using a common form of the Mori-Tanaka equation.
    # This form shows the effective modulus as the matrix modulus plus a correction term.
    # We use Vm (matrix volume fraction) in the inverse term for clarity,
    # remembering the physical constraint that Vm + Vf = 1.
    # The term `(Vm * I + Vf * A)**-1` represents the inverse of the tensor in parentheses.
    inverse_term = (Vm * I + Vf * A)**-1
    correction_term = Vf * (Cf - Cm) * A * inverse_term
    expression = Cm + correction_term

    # Print the final equation.
    # The request to "output each number in the final equation" for a symbolic
    # problem is interpreted as printing the complete symbolic equation clearly,
    # showing how C is constructed from all the given quantities.
    print("The expression for the effective average elastic moduli C in the Mori-Tanaka model is:")

    # We create a sympy Equation object for pretty printing.
    final_equation = sympy.Eq(C, expression)

    # Using sympy.pprint for a more readable, formatted output.
    # If this script is run in a terminal that doesn't support unicode,
    # use_unicode=False is a safe choice.
    sympy.pprint(final_equation, use_unicode=False)

if __name__ == '__main__':
    mori_tanaka_expression()