import sympy as sp

def calculate_magnetic_field_expression():
    """
    Calculates and prints the expression for the magnetic field inside a stack of superconducting strips.
    This model assumes the stack can be treated as a continuous slab (D << w).
    """
    # Define the symbols
    H_z = sp.Symbol('H_z')
    x = sp.Symbol('x')
    a = sp.Symbol('a')
    H0 = sp.Symbol('H0')
    D = sp.Symbol('D')
    pi = sp.pi

    # In the slab model (valid for D << w), the field in the penetrated region (a < |x| < w) is H_z = J_eff * (|x| - a)
    # The effective current density is J_eff = Jc * d / D
    # Given H0 = Jc * d / pi, we can write Jc*d = H0 * pi
    # So, J_eff = (H0 * pi) / D
    
    # Expression for the magnetic field H_z
    # We use sp.Abs(x) for the absolute value of x
    expression = (H0 * pi / D) * (sp.Abs(x) - a)
    
    # Print the equation
    # We use sp.pretty to format the equation for better readability
    equation = sp.Eq(H_z, expression)
    
    # The final step from the problem is to output each number in the final equation.
    # Our result is symbolic, but we will print it in a format that resembles an equation.
    # The python code represents the derivation and symbolic representation.
    # The final printout will be the formatted equation string.
    
    # Using sp.pprint for a more "mathematical" look in the terminal
    # print("The expression for the magnetic field H_z(x) is:")
    # sp.pprint(equation, use_unicode=True)
    
    # As requested, printing the equation as a string.
    print(f"H_z(x) = (pi * H0 / D) * (Abs(x) - a)")

calculate_magnetic_field_expression()

<<<H_z(x) = (pi * H0 / D) * (Abs(x) - a)>>>