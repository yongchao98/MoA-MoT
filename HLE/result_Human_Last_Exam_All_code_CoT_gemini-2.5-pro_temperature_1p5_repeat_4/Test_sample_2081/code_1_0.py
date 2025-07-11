import math

def solve_for_R_and_print_equation():
    """
    This function calculates the radius R based on the solvability condition
    of the given boundary-value problem and prints the results.
    """
    
    # Step 1: Define the given parameter T and calculate e^T.
    # T = ln(10^34), so e^T = 10^34.
    # We use this direct value to maintain precision.
    eT = 10.0**34
    
    # Step 2: Calculate R² using the simplified formula derived from the
    # solvability condition.
    # R² = 0.5 * (e^T + 1) * e^T
    R_squared = 0.5 * (eT + 1) * eT
    
    # Step 3: Calculate R, the radius of the sphere of valid initial conditions.
    R = math.sqrt(R_squared)
    
    # Step 4: Print the final equation and its numerical components as requested.
    # The equation defines the set of initial values (x₀, y₀, z₀) for which solutions exist.
    print("The solvability condition results in the following equation for the initial values (x₀, y₀, z₀):")
    # We print the numbers 1, 1, 1, and the calculated value for R².
    print(f"1 * x₀² + 1 * y₀² + 1 * z₀² = {R_squared}")
    
    # Step 5: Print the final answer for R.
    print(f"\nThe radius R of the sphere of valid initial values is:")
    print(R)

# Execute the function to get the solution.
solve_for_R_and_print_equation()