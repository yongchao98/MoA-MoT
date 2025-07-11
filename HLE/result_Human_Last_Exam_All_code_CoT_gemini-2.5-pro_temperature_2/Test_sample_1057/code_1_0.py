import math

def solve_joule_heat():
    """
    This function prints the derived formula for the total Joule heat dissipated.
    The problem is symbolic, so the output is the final equation rather than a number.
    """
    
    # Define symbolic representations of the variables from the problem
    initial_radius = "a"
    initial_potential = "V"
    
    # The physical constants involved
    vacuum_permittivity = "ε₀"
    pi_symbol = "π"
    
    # The derivation leads to the formula H = 4 * pi * epsilon_0 * a * V^2.
    # We will print this formula, highlighting its components as per the instructions.

    print("The problem asks for the total Joule heat (H) dissipated during the process.")
    print("Based on the principles of energy conservation, we derived the final formula.")
    print("\nThe final equation for the total Joule heat dissipated is:")
    
    # We print each part of the formula clearly.
    # The prompt asks to "output each number in the final equation", which is '4' and '2' (for the exponent).
    final_equation = f"H = 4 * {pi_symbol} * {vacuum_permittivity} * {initial_radius} * {initial_potential}²"
    print(final_equation)
    
    print("\nWhere the components of the equation are:")
    print("  - The number '4' is a numerical coefficient from the derivation.")
    print(f"  - '{pi_symbol}' is the mathematical constant pi (approx. {math.pi}).")
    print(f"  - '{vacuum_permittivity}' is the physical constant for the vacuum permittivity.")
    print(f"  - '{initial_radius}' is the initial radius of the sphere.")
    print(f"  - '{initial_potential}²' is the initial potential of the sphere, squared.")

solve_joule_heat()