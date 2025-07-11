import sympy as sp
from sympy import sin, cosh, sinh

def solve_force_equation():
    """
    This function uses the sympy library to symbolically represent the final
    derived equation for the force per unit area and prints its components.
    """
    # Define the symbols used in the equation
    K0, a, y, d, mu0, mu = sp.symbols('K_0 a y d mu_0 mu', real=True, positive=True)

    # Define a symbol for the unit vector in the x-direction
    i_x_hat = sp.Symbol('\\hat{i}_x')

    # The problem asks for the force per unit area on the x=d interface.
    # The derived formula for this force is:
    # f/area = (mu_0 / 2) * (K_0**2 * sin(a*y)**2) / (cosh(a*d) + (mu_0/mu)*sinh(a*d))**2 * i_x_hat

    # Let's build this expression step by step, as requested.
    # "each number in the final equation" is interpreted as each component part.

    print("The final equation for the force per unit area is constructed from the following components:\n")

    # Component 1: The coefficient
    coefficient = mu0 / 2
    print(f"Coefficient part: {coefficient}")

    # Component 2: The source current term in the numerator
    numerator_current_term = K0**2
    print(f"Numerator current term: {numerator_current_term}")

    # Component 3: The spatial variation term in the numerator
    numerator_spatial_term = sin(a*y)**2
    print(f"Numerator spatial term: {numerator_spatial_term}")

    # The full numerator is the product of these parts
    full_numerator = coefficient * numerator_current_term * numerator_spatial_term
    print(f"--- Full Numerator: {full_numerator}")


    # Component 4: The first term inside the denominator's parenthesis
    denominator_term_1 = cosh(a*d)
    print(f"\nDenominator component 1 (from geometry): {denominator_term_1}")

    # Component 5: The second term inside the denominator's parenthesis
    denominator_term_2 = (mu0/mu) * sinh(a*d)
    print(f"Denominator component 2 (from materials and geometry): {denominator_term_2}")

    # The full denominator is the square of the sum of these parts
    full_denominator = (denominator_term_1 + denominator_term_2)**2
    print(f"--- Full Denominator: {full_denominator}")

    # Component 6: The direction vector
    direction_vector = i_x_hat
    print(f"\nDirection: {direction_vector}")


    # Now, let's assemble and print the final, complete formula.
    force_per_area_expression = (full_numerator / full_denominator) * direction_vector

    print("\n--------------------------------------------------")
    print("Final Expression for Force Per Unit Area:")
    print("--------------------------------------------------")
    # Using sympy's pretty print for a more readable output
    sp.pprint(force_per_area_expression, use_unicode=True)
    print("--------------------------------------------------")

if __name__ == '__main__':
    solve_force_equation()