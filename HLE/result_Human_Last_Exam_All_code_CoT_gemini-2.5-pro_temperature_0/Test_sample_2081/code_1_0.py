import math
from decimal import Decimal, getcontext

def solve_for_R():
    """
    Calculates the radius R of the sphere of initial conditions for which
    the nonlinear boundary value problem has a solution.
    """
    # Set a high precision for the calculation to handle large numbers accurately.
    getcontext().prec = 100

    # From the problem statement, T = ln(10^34), which implies e^T = 10^34.
    # We use the Decimal type for high-precision arithmetic.
    e_T = Decimal(10)**34

    # The solvability condition for the system leads to the equation for the 
    # set of initial values (x₀, y₀, z₀), which forms a sphere:
    # 1*x₀² + 1*y₀² + 1*z₀² = R²
    #
    # The squared radius R² is derived from the solvability condition as:
    # R² = (1/2) * e^T * (e^T + 1)
    
    # Define the coefficients and constants in the calculation
    coeff_x_sq = 1
    coeff_y_sq = 1
    coeff_z_sq = 1
    half = Decimal('0.5')

    # Calculate R²
    R_squared = half * e_T * (e_T + 1)

    # Calculate R by taking the square root
    R = R_squared.sqrt()

    # As requested, we output each number in the final equation for the sphere.
    print("The final equation for the set of initial values is a sphere: c_x*x₀² + c_y*y₀² + c_z*z₀² = R²")
    print(f"Coefficient c_x: {coeff_x_sq}")
    print(f"Coefficient c_y: {coeff_y_sq}")
    print(f"Coefficient c_z: {coeff_z_sq}")
    print(f"The value of R² is: {R_squared}")
    print(f"The radius R is: {R}")

solve_for_R()