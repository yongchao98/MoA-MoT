import sympy as sp

def solve_sliding_rod():
    """
    This function derives and prints the expression for the angle theta
    at which the rod begins to slide.
    """
    # Define the symbolic variables for the problem
    L, l, mu, theta = sp.symbols('L l mu theta', real=True, positive=True)

    # Numerator of the expression for tan(theta) / mu
    numerator = L**2 / 12 + 2 * l**2
    
    # Denominator of the expression for tan(theta) / mu
    denominator = L**2 / 12 + 3 * l**2
    
    # Expression for tan(theta)
    tan_theta_expr = mu * (numerator / denominator)
    
    # The equation for theta
    # We can write this as an explicit equation for theta using atan
    theta_solution = sp.atan(tan_theta_expr)
    
    # Print the final result in a readable format
    # The printout shows all numbers (12, 2, 3) in the final equation.
    print("The angle theta at which the rod begins to slide is given by the expression:")
    print(f"theta = {theta_solution}")

# Run the function to display the result
solve_sliding_rod()