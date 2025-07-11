import numpy as np
from scipy.integrate import quad

def main():
    """
    This script calculates the value of the time-averaged integral from the problem.
    The derivation shows that the problem reduces to calculating the definite integral:
    I = integral from 0 to infinity of (4 / (3 * exp(2*t^2) + 1)) dt.
    We use numerical integration to find the value of I.
    """

    # Define the integrand function
    # The variable of integration is tau (t in the code)
    integrand = lambda t: 4 / (3 * np.exp(2 * t**2) + 1)

    # Perform the numerical integration from 0 to infinity
    integral_value, error = quad(integrand, 0, np.inf)

    # The final answer is the value of the integral
    final_answer = integral_value

    # The problem asks to output each number in the final equation.
    # The final equation is I = integral_value.
    # Let's print the components for clarity, although it's straightforward here.
    
    print("The integral to be calculated is: I = integral_from_0_to_inf[4 / (3*exp(2*t^2) + 1)] dt")
    print(f"The numerical value of the integral is: {final_answer:.10f}")
    print(f"The estimated error of the calculation is: {error:.2e}")
    
    # Final answer as requested by the format
    print("\nFinal Answer:")
    print(final_answer)


if __name__ == "__main__":
    main()
