import numpy as np
from scipy import integrate

def evaluate_the_integral():
    """
    This function evaluates the definite integral numerically.
    The integral is I = ∫[from 0 to φ³-1] Re[1 / (1 + e^(arctan(ln(cos(x/e)))))^i] dx.
    """
    
    # Step 1: Define the necessary mathematical constants.
    # φ (phi) is the golden ratio.
    phi = (1 + np.sqrt(5)) / 2
    # e is Euler's number.
    e = np.e

    # Step 2: Define the limits of integration.
    # The lower limit 'a' is 0.
    lower_limit = 0.0
    # The upper limit 'b' is φ³ - 1, which simplifies to 2*φ.
    upper_limit = phi**3 - 1

    # Step 3: Define the integrand function.
    # As derived in the plan, the complex expression simplifies to a real-valued function.
    # f(x) = cos(ln(w(x))), where w(x) = 1 + exp(arctan(ln(cos(x/e)))).
    # This simplified form is more efficient and numerically stable.
    def integrand(x):
        # For the given integration range [0, φ³-1], the value x/e is always
        # within [0, π/2), ensuring cos(x/e) is positive and its logarithm is defined.
        cos_val = np.cos(x / e)
        log_val = np.log(cos_val)
        arctan_val = np.arctan(log_val)
        exp_val = np.exp(arctan_val)
        w = 1 + exp_val
        
        return np.cos(np.log(w))

    # Step 4: Perform the numerical integration using scipy.integrate.quad.
    # The function returns the result of the integral and an estimate of the error.
    result, error_estimate = integrate.quad(integrand, lower_limit, upper_limit)

    # Step 5: Print the final results as requested.
    # The final equation is: ∫[from lower_limit to upper_limit] f(x) dx = result.
    # The instruction is to output each number in this equation.
    print("The final equation is of the form: Integral = ∫ f(x) dx from a to b")
    print("\nThe numbers in the final equation are:")
    print(f"Lower limit a = {lower_limit}")
    print(f"Upper limit b = {upper_limit}")
    print(f"Result of the integral = {result}")
    print(f"(Estimated error of the calculation is {error_estimate})")

# Run the evaluation function
evaluate_the_integral()