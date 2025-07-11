import math

def solve_for_R():
    """
    This function calculates the radius R based on the derived solvability condition.
    The set of initial values (x_0, y_0, z_0) for which solutions to the
    nonlinear problem exist forms a sphere defined by x_0^2 + y_0^2 + z_0^2 = R^2.
    """
    
    # Given T = ln(10^34), we can find e^T directly.
    # Using floating-point numbers for the calculation.
    e_T = 10.0**34
    
    # The radius squared R^2 is given by the formula:
    # R^2 = 0.5 * e^T * (e^T + 1)
    # The +1 is negligible in standard float precision but we include it.
    R_squared = 0.5 * e_T * (e_T + 1.0)
    
    # The radius R is the square root of R^2.
    R = math.sqrt(R_squared)
    
    # Print the final equation for the initial values
    print("The solvability condition imposes a constraint on the initial values (x_0, y_0, z_0).")
    print("This constraint describes a sphere with the following equation:")
    print(f"x_0^2 + y_0^2 + z_0^2 = R^2")
    print(f"where R^2 = {R_squared:.4e}")
    print("\nExecuting the calculation...")
    
    # Output each number in the final equation.
    # The coefficients of x_0^2, y_0^2, z_0^2 are 1.
    print(f"Coefficient of x_0^2: 1")
    print(f"Coefficient of y_0^2: 1")
    print(f"Coefficient of z_0^2: 1")
    print(f"Value of R^2: {R_squared}")
    
    # Print the final result for R.
    print(f"\nThe radius of the sphere is R = {R}")

solve_for_R()