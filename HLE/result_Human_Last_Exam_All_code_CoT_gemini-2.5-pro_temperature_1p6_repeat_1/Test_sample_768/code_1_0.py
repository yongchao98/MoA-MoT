import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    Evaluates the definite integral using numerical methods.
    """
    # Define the golden ratio, phi
    phi = (1 + np.sqrt(5)) / 2

    # Set the integration limits
    lower_limit = 0
    upper_limit = phi**3 - 1

    # Define the integrand function f(x)
    # f(x) = Re[1 / (1 + exp(arctan(ln(cos(x/e)))))^i]
    # which simplifies to cos(ln(1 + exp(arctan(ln(cos(x/e))))))
    def integrand(x):
        # The argument of ln, cos(x/e), is positive over the integration interval
        cos_val = np.cos(x / np.e)
        
        # Calculate A(x) = 1 + exp(...)
        A = 1 + np.exp(np.arctan(np.log(cos_val)))
        
        # The integrand is cos(ln(A))
        return np.cos(np.log(A))

    # Perform the numerical integration
    # The quad function returns the result and an estimated error
    result, error = quad(integrand, lower_limit, upper_limit, limit=200)

    # The prompt requests to output the result using the print function.
    # The numerical result is expected to be very close to phi**2.
    # phi**2 = phi + 1 = (3 + sqrt(5))/2 approx 2.6180339887
    # "Remember in the final code you still need to output each number in the final equation!"
    # The final equation is: Integral = result. We print the numbers involved.
    print(f"The integral from {lower_limit} to {upper_limit:.10f} is: {result:.10f}")
    # The problem asks for evaluation. The numerical result is the evaluation.
    # We can also compare it to the theoretical value phi^2.
    print(f"The value of phi^2 is: {phi**2:.10f}")

solve_integral()