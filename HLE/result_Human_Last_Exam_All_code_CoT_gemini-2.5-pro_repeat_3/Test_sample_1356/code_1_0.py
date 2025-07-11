import sympy

def solve_pendulum_period():
    """
    This function calculates and prints the symbolic period of the described physical system.
    
    The final formula for the period T is derived using Lagrangian mechanics for small oscillations.
    It is expressed in terms of the disk radius R and the acceleration due to gravity g.
    
    The formula is: T = 2 * pi * sqrt( (101 * R) / (40 * g) )
    """
    
    # Define the symbolic variables for the equation
    R, g, pi = sympy.symbols('R g pi')
    
    # The numerical coefficients derived from the Lagrangian mechanics analysis
    numerator_coefficient = 101
    denominator_coefficient = 40
    two = 2
    
    # Build the symbolic expression for the period T
    period_expression = two * pi * sympy.sqrt((numerator_coefficient * R) / (denominator_coefficient * g))
    
    # Print the final equation clearly, showing each number
    print("The period of motion, T, is given by the following equation:")
    
    # We use the integer values directly in the printout as requested
    print(f"T = {two} * π * sqrt( ({numerator_coefficient} * R) / ({denominator_coefficient} * g) )")

# Execute the function to get the answer
solve_pendulum_period()