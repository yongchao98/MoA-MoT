import sympy

def calculate_total_mass():
    """
    This function symbolically calculates the final answer based on the problem's parameters.
    """
    # Define the symbolic variables from the problem statement.
    # q is the size of the constant field of the global function field.
    # q_v is the order of the residual field of the local field.
    q = sympy.Symbol('q')
    q_v = sympy.Symbol('q_v')

    # Based on the analysis, the mass of the measure mu is 1 / (q - 1).
    mass_mu = 1 / (q - 1)

    # This is the coefficient given in the problem.
    coefficient = (q_v * (q - 1)) / (q_v - 1)

    # The final quantity is the product of the coefficient and the mass.
    total_mass = coefficient * mass_mu

    # Use sympy to simplify the expression.
    simplified_mass = sympy.simplify(total_mass)

    # The simplified result is q_v / (q_v - 1).
    # We will now print the final equation as requested, showing each component.
    numerator, denominator = simplified_mass.as_numer_denom()
    
    # The denominator is an expression (q_v - 1). We extract its parts.
    den_variable = denominator.args[0]
    den_constant = -denominator.args[1]

    print(f"Total Mass = {numerator} / ({den_variable} - {den_constant})")

calculate_total_mass()