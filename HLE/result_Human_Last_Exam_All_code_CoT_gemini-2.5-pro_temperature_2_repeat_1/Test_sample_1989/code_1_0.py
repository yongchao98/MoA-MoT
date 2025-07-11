import sympy as sp

def display_corrector_formula():
    """
    This function defines and displays the corrector factor for the given PDE problem using SymPy.
    """
    # Define the symbolic variables
    r, theta, A, B = sp.symbols('r, theta, A, B', real=True, positive=True)

    # The original asymptotic behavior without the corrector
    # base_omega = (1 / sp.sqrt(r)) * sp.exp(-r * (1 - sp.cos(theta)))

    # Define the exponent for the 'r' term in the corrector
    power_exponent = A * (1 - sp.cos(theta)) + B * sp.sin(theta)

    # Define the terms of the corrector
    term1_power_of_r = r**power_exponent
    term2_exp_of_theta = sp.exp(B * theta)

    # The full corrector factor
    corrector = term1_power_of_r * term2_exp_of_theta
    
    print("The corrector to the large-distance behavior of omega is:")
    
    # We use SymPy's pretty print for a nice mathematical layout
    sp.init_printing(use_unicode=True)
    sp.pprint(corrector)

    # To satisfy the requirement of outputting each number in the equation,
    # let's break down the expression and show its components.
    print("\n--- Breakdown of the Corrector's expression ---")
    print("\nThe corrector is a product of two terms: C = (Term 1) * (Term 2)\n")

    print("Term 1: A power of r")
    sp.pprint(term1_power_of_r)
    print("The exponent of r is:")
    sp.pprint(power_exponent)
    print("This exponent contains the number 1.")


    print("\nTerm 2: An exponential function of theta")
    sp.pprint(term2_exp_of_theta)
    print("The exponent is:")
    sp.pprint(B * theta)
    

display_corrector_formula()