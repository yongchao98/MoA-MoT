import math

def solve_for_R():
    """
    This function calculates the radius R based on the solvability condition
    of the given boundary-value problem.
    """
    
    # From the problem, T = ln(10^34), so exp(T) = 10^34.
    # We use floating-point numbers for calculations involving large magnitudes.
    e_T = 10.0**34
    
    # The solvability condition for the nonlinear problem is that the set of initial
    # values (x_0, y_0, z_0) must lie on a sphere x_0^2 + y_0^2 + z_0^2 = R^2.
    # The radius squared R^2 is given by the formula:
    # R^2 = alpha / (1 - exp(-T))
    # where alpha = (exp(2*T) - 1) / 2.

    # We can calculate the components of this formula first.
    # Calculate alpha:
    e_2T = e_T**2
    alpha = (e_2T - 1) / 2
    
    # Calculate the denominator:
    e_minus_T = 1.0 / e_T
    denominator = 1 - e_minus_T

    # To improve numerical stability, we can simplify the expression for R^2:
    # R^2 = [0.5 * (e^T - 1)(e^T + 1)] / [(e^T - 1) / e^T]
    # R^2 = 0.5 * e^T * (e^T + 1)
    # We will use this stable formula for the final calculation.
    R_squared = 0.5 * e_T * (e_T + 1)
    R = math.sqrt(R_squared)

    # Now, we print the results as requested, showing each number in the final equation.
    print("The set of initial values (x_0, y_0, z_0) for which solutions exist is a sphere defined by the equation:")
    print("x_0^2 + y_0^2 + z_0^2 = R^2")
    print("\nTo find R, we use the formula R^2 = alpha / (1 - exp(-T))")
    print("\nFirst, we calculate the values for alpha and the denominator:")
    print(f"alpha = (exp(2*T) - 1) / 2 = ({e_2T:.4e} - 1) / 2 = {alpha:.4e}")
    print(f"1 - exp(-T) = 1 - {e_minus_T:.4e} = {denominator}")
    
    print("\nSubstituting these values, the final equation for the initial values is:")
    print(f"x_0^2 + y_0^2 + z_0^2 = {alpha:.4e} / {denominator}")
    print(f"x_0^2 + y_0^2 + z_0^2 = {R_squared:.4e}")

    print("\nThe radius R is the square root of this value:")
    print(f"R = {R:.4e}")

if __name__ == "__main__":
    solve_for_R()
