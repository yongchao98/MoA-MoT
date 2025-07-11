import sympy

def generate_mori_tanaka_expression():
    """
    This function uses the sympy library to generate and print the symbolic expression
    for the effective elastic moduli C of a composite material based on the Mori-Tanaka model.
    """
    # Define scalar symbols for volume fractions. Vm + Vf = 1.
    Vf = sympy.Symbol('Vf')
    Vm = sympy.Symbol('Vm')

    # Define non-commuting symbols for the fourth-order tensors.
    # The order of multiplication matters for tensors, hence commutative=False.
    I = sympy.Symbol('I', commutative=False)
    Cf = sympy.Symbol('Cf', commutative=False)
    Cm = sympy.Symbol('Cm', commutative=False)
    A = sympy.Symbol('A', commutative=False)
    C = sympy.Symbol('C', commutative=False)

    # Build the expression for C using the second derived form, which is common in literature.
    # C = Cm + Vf * (Cf - Cm) * A * inverse(Vm*I + Vf*A)
    
    # Term representing the difference in stiffness between fiber and matrix
    stiffness_difference = Cf - Cm
    
    # Term inside the inverse
    inverse_term = Vm * I + Vf * A
    
    # The full expression for the effective stiffness tensor C
    expression_for_C = Cm + Vf * stiffness_difference * A * inverse_term**-1

    # Print the final equation. This fulfills the requirement to show each term.
    # The sympy output is a standard representation of the formula.
    print(f"The expression for the effective elastic moduli C is:")
    print(f"{C} = {expression_for_C}")

if __name__ == '__main__':
    generate_mori_tanaka_expression()