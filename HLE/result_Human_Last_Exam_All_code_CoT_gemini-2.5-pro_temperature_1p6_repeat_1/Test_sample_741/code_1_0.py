import numpy as np
from scipy.special import iv
from scipy.optimize import root_scalar

# The problem is to find the largest x for which the summation converges to 0.
# The summation can be expressed in terms of the modified Bessel function of the first kind:
# sum = I_{x-1}(2)
# We need to find the largest root x for the equation I_{x-1}(2) = 0.
# Let nu = x - 1. We are looking for the largest root of I_nu(2) = 0.

# Define the function f(nu) = I_nu(2), whose root we want to find.
def bessel_function_of_order(nu):
    """
    Computes the modified Bessel function of the first kind, I_nu(2).
    """
    return iv(nu, 2.0)

# From the properties of Bessel functions, we know the roots for the order 'nu'
# are real and negative. We seek the largest (least negative) root.
# A plot or prior knowledge suggests this root is in the interval [-3, -2].
# We use the root_scalar function to find the root within this bracket.
try:
    solution = root_scalar(bessel_function_of_order, bracket=[-3, -2])
    nu_root = solution.root
    
    # Calculate x from the relationship nu = x - 1
    x_value = nu_root + 1
    
    # Print the final x value in the specified format {-a.bbb}
    # The final equation is x = value, and we print the number for this equation.
    print(f"{{{x_value:.3f}}}")

except (ImportError, ValueError) as e:
    print("An error occurred. Please ensure you have scipy installed (`pip install scipy`).")
    print(f"Error details: {e}")
