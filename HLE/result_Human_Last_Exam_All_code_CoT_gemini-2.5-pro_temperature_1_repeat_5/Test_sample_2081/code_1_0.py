import math

def solve_for_R():
    """
    This function calculates the radius R of the sphere of initial conditions
    for the given boundary-value problem.
    """
    # Step 1: Define the given constants from the problem.
    # We are given T = ln(10^34), which implies e^T = 10^34.
    # The problem specifies alpha in terms of T, but the solvability condition
    # simplifies to a direct relationship between R and T.
    e_T = 10.0**34

    # Step 2: Use the derived formula for the square of the radius R.
    # The solvability condition for the system leads to the equation for the
    # set of initial values (x0, y0, z0):
    # x0^2 + y0^2 + z0^2 = R^2
    # where R^2 = alpha / (1 - e^(-T)).
    # Substituting alpha = 0.5 * (e^(2T) - 1) and simplifying gives:
    # R^2 = 0.5 * (e^T + 1) * e^T
    
    # The numbers in this final equation are:
    c1 = 0.5
    c2 = 1.0

    print("The derived equation for the square of the radius R is:")
    print(f"R^2 = {c1} * (e^T + {c2}) * e^T")
    print("-" * 30)

    # Step 3: Substitute the value of e^T and calculate R^2.
    # Note: For very large e_T, (e_T + 1) is computationally equal to e_T
    # in standard floating-point arithmetic.
    R_squared = c1 * (e_T + c2) * e_T

    print(f"Substituting e^T = {e_T:.0e}:")
    # We show the components of the calculation as requested.
    print(f"R^2 = {c1} * ({e_T:.0e} + {c2}) * {e_T:.0e}")
    print(f"R^2 = {R_squared:.4e}")
    print("-" * 30)

    # Step 4: Calculate R by taking the square root.
    R = math.sqrt(R_squared)

    # Step 5: Print the final result for R.
    print("The value of R is the square root of R^2:")
    print(f"R = {R:.15e}")

solve_for_R()