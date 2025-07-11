import numpy as np
from scipy.integrate import quad

def solve_temperature():
    """
    This function calculates the temperature T at time t = arccosh(2)/2
    by numerically evaluating the analytical solution derived from the
    differential equation.
    The final expression for the temperature T at a transformed variable u=cosh(2t) is:
    T(u) = -u + exp(u**2) / (exp(1) + integral from 1 to u of exp(s**2) ds)
    We evaluate this at u = 2.
    """

    # Define the integrand for the numerical integration
    def integrand(s):
        return np.exp(s**2)

    # The problem requires evaluation at a specific time t, which corresponds to u=2
    u_final = 2.0

    # These are the components of the final equation to calculate the temperature
    term_A = -u_final
    numerator_val = np.exp(u_final**2)
    denominator_addend_val = np.exp(1)

    # Perform the numerical integration for the integral part of the denominator
    # The integral is evaluated from u=1 (at t=0) to u=2 (at the final t)
    integral_val, _ = quad(integrand, 1, u_final)

    # Calculate the full denominator
    denominator_val = denominator_addend_val + integral_val

    # Calculate the final temperature
    temperature = term_A + numerator_val / denominator_val
    
    # As requested, output the numbers in the final equation
    print("The final temperature is calculated using the formula: T = A + B / (C + I)")
    print(f"Value of A (which is -u): {term_A}")
    print(f"Value of B (which is exp(u^2)): {numerator_val}")
    print(f"Value of C (which is exp(1)): {denominator_addend_val}")
    print(f"Value of the integral I from 1 to 2: {integral_val}")
    print("\nFinal calculation:")
    print(f"T = {term_A} + {numerator_val} / ({denominator_addend_val} + {integral_val})")
    print(f"T = {temperature}")

solve_temperature()

# The final answer in the required format
# We'll calculate the value again to be sure
u = 2.0
term_A = -u
numerator = np.exp(u**2)
denominator_addend = np.exp(1)
integral, err = quad(lambda s: np.exp(s**2), 1, u)
denominator = denominator_addend + integral
final_temp = term_A + numerator / denominator
print(f"<<<{final_temp:.3f}>>>")