import sympy

def solve_steady_state_pi0():
    """
    This function derives and prints the steady-state probability pi_0
    for the given birth-death process in terms of rho.
    The derivation leads to the equation pi_0 = exp(-rho).
    This script will print the final equation symbolically.
    """
    
    # Define pi_0 and rho as symbolic variables for clear representation
    pi_0 = sympy.Symbol('pi_0')
    rho = sympy.Symbol('rho')
    
    # The derived expression for pi_0 is exp(-rho).
    # We can create a symbolic equation to represent this relationship.
    final_equation = sympy.Eq(pi_0, sympy.exp(-rho))
    
    # Print the final derived equation.
    # The equation itself contains all the elements (variables, constants, operators).
    print("The final equation for the steady-state probability pi_0 in terms of rho is:")
    print(final_equation)

solve_steady_state_pi0()