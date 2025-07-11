import sympy

def solve_force_equation():
    """
    This function symbolically represents the formula for the x-directed total force
    on the specified conducting material.
    """
    # Define the symbols used in the equation
    a, D, mu_0, I_0, sigma_1, sigma_2 = sympy.symbols('a D mu_0 I_0 sigma_1 sigma_2')
    
    # The problem asks for the force on the material in the region s < x < 2s (block 2).
    # Based on analyzing the answer choices for physical consistency, Option A is the most plausible.
    # The formula for the force is:
    # F_x = -a*D * (mu_0/2) * (I_0^2/D^2) * (sigma_2/(sigma_1 + sigma_2))^2
    
    # Let's break down the expression to show the components, including the numeric parts.
    
    # The coefficient is -1/2
    coefficient = sympy.Rational(-1, 2)
    
    # The geometric and physical constant terms
    term_constants = mu_0 * a * D
    
    # The source current term
    term_current = (I_0 / D)**2
    
    # The conductivity-dependent term (related to current division)
    term_conductivity = (sigma_2 / (sigma_1 + sigma_2))**2

    # Combine the terms to form the final expression for the force
    F_x = coefficient * term_constants * term_current * term_conductivity
    
    # Print the final equation clearly showing each component
    print("The x-directed total force on the conducting material is given by the equation:")
    
    # The instruction asks to output each number. We will print the components of the formula.
    # Coefficient part
    print("Coefficient = -1/2")
    
    # Geometric and permeability part
    print("Geometric and constants part = a * D * mu_0")
    
    # Current part
    print("Current term = (I_0/D)^2")
    
    # Conductivity part
    print("Conductivity term = (sigma_2 / (sigma_1 + sigma_2))^2")

    # Full equation reconstruction
    print("\nFinal Equation:")
    # Using sympy.pretty_print for a clearer mathematical output
    final_eq_str = f"F_x = {sympy.pretty(F_x, use_unicode=False)}"
    print(final_eq_str.replace('mu_0', 'μ₀').replace('sigma_1', 'σ₁').replace('sigma_2', 'σ₂'))


solve_force_equation()