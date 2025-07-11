import sympy
from sympy import symbols, Eq, Mul, Pow, Rational

def solve_force_equation():
    """
    This function prints the derived formula for the x-directed total force.
    """
    # Define the symbolic variables for the system parameters
    F_x = symbols('F_x')
    a, D, mu_0, I_0, sigma_1, sigma_2 = symbols('a D mu_0 I_0 sigma_1 sigma_2')

    # Construct the terms of the equation based on the derivation
    # Coefficient term
    coeff = Mul(Rational(-1, 2), a, D, mu_0, evaluate=False)
    
    # Current term
    current_term = Pow(I_0, 2) / Pow(D, 2)
    
    # Conductivity term
    conductivity_term = Pow(sigma_2 / (sigma_1 + sigma_2), 2)
    
    # Full expression for the force
    rhs = coeff * current_term * conductivity_term
    
    # Create the equation
    equation = Eq(F_x, rhs)

    # Print the final equation in a readable format
    print("The derived x-directed total force on the conducting material is:")
    
    # We use sympy's pretty print for a clear mathematical representation
    sympy.init_printing(use_unicode=True)
    print(equation)
    
    print("\nBreaking down the final equation from choice A:")
    print(f"F_x = (Term 1) * (Term 2) * (Term 3)")
    print(f"Term 1 (Coefficient): -a*D * mu_0 / 2")
    print(f"Term 2 (Current part): (I_0^2 / D^2)")
    print(f"Term 3 (Conductivity part): (sigma_2 / (sigma_1 + sigma_2))^2")
    
    final_expression_str = f"F_x = -a*D * (mu_0/2) * (I_0^2/D^2) * (sigma_2/(sigma_1 + sigma_2))^2"
    print(f"\nFinal Equation String: {final_expression_str}")

solve_force_equation()