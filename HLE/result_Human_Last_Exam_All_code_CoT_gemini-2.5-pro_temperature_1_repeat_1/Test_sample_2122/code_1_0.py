import numpy as np
from scipy.integrate import quad

def solve_temperature():
    """
    Solves for the temperature T at the specified time by numerically evaluating
    the integral in the solution of the differential equation.
    """
    # The integrand in the solution for T(x)
    def integrand(s):
        return np.exp(s**2)

    # The problem is transformed from variable t to x.
    # The final time corresponds to x = 2.
    x_final = 2
    
    # We need to compute the integral from x=1 to x=2
    integral_start = 1
    integral_end = 2
    
    # Numerically calculate the definite integral
    # I = integral from 1 to 2 of exp(s^2) ds
    integral_value, _ = quad(integrand, integral_start, integral_end)
    
    # Constants needed for the final calculation
    e_val = np.exp(1)
    e_pow_4_val = np.exp(x_final**2)
    
    # The solution for T(x=2) is T(2) = -2 + exp(4) / (exp(1) + I)
    numerator = e_pow_4_val
    denominator = e_val + integral_value
    
    # Calculate the final temperature
    final_temperature = -x_final + numerator / denominator

    # Output the components and the final result
    print("The final temperature is found by evaluating the expression:")
    print(f"T(x={x_final}) = -{x_final} + exp({x_final}^2) / (exp(1) + I)")
    print(f"where I is the integral of exp(s^2) from {integral_start} to {integral_end}.\n")
    
    print("The values used in the final equation are:")
    print(f"Value of the integral I = {integral_value}")
    print(f"Value of exp(1) = {e_val}")
    print(f"Value of exp({x_final}^2) = {e_pow_4_val}\n")

    print("Substituting these numbers into the equation:")
    print(f"T(2) = -{x_final} + {numerator} / ({e_val} + {integral_value})")
    print(f"T(2) = -{x_final} + {numerator / denominator}")
    print(f"The temperature at the specified time is T = {final_temperature}")

solve_temperature()