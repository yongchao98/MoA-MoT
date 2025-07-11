import sympy as sp

def print_force_formula():
    """
    This function prints the symbolic formula for the force per unit area on the x=d interface.
    The formula corresponds to one of the multiple-choice answers.
    """
    mu_0, mu, K_0, a, y, d = sp.symbols('mu_0 mu K_0 a y d')
    
    # Define the numerator of the force expression
    numerator = mu_0 / 2 * (K_0**2 * sp.sin(a*y)**2)
    
    # Define the denominator of the force expression
    denominator_base = sp.cosh(a*d) + (mu_0 / mu) * sp.sinh(a*d)
    denominator = denominator_base**2
    
    # The force is a vector in the x-direction
    force_direction = "\\hat{i}_x"
    
    # We will print the formula in a readable format.
    # Using sp.pretty_print for a nicer mathematical layout.
    
    print("The force per unit area is given by:")
    
    # To construct the pretty print output, we can use a string and replace parts
    # since SymPy's pretty printer for fractions can be verbose.
    
    force_str = f"""
         mu_0 * K_0^2 * sin^2(a*y)
f/area = -------------------------------- * {force_direction}
            2 * [cosh(a*d) + (mu_0/mu)*sinh(a*d)]^2
"""
    
    # Print the final expression
    # Due to the complexity of printing a nicely formatted equation as a string,
    # let's construct it term by term for clarity.
    
    print(f"f/area = (mu_0 / 2) * (K_0**2 * sin(a*y)**2) / ( (cosh(a*d) + (mu_0/mu)*sinh(a*d))**2 )  in the x-direction")

    # The problem asks to output the formula itself.
    # The prompt says: "Remember in the final code you still need to output each number in the final equation!"
    # The formula is symbolic, so we will print the variable names.

    final_expression_latex = r"\frac{\vec{f}}{\text{area}} = \frac{\mu_0}{2} \frac{K_0^2 \sin^2(ay)}{\left[ \cosh(ad) + \frac{\mu_0}{\mu} \sinh(ad) \right]^2} \hat{i}_x"
    
    print("\nFinal formula:")
    print(final_expression_latex)


print_force_formula()