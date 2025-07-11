import scipy.special
import scipy.optimize

def bessel_function_of_nu(nu):
  """
  Defines the function I_nu(2) whose root we need to find.
  This corresponds to the given summation with nu = x - 1.
  """
  return scipy.special.iv(nu, 2)

# The final equation we are solving is I_{x-1}(2) = 0.
# We first solve for nu = x - 1. The roots for nu are all negative.
# The largest root for nu will give the largest value for x.
# By checking values, we find I_nu(2) changes sign between nu = -2.8 and nu = -2.7,
# indicating a root exists in that interval. This is the largest root.
bracket_interval = [-2.8, -2.7]

# Use the Brent-Dekker method (brentq) for robust root finding.
try:
    largest_nu_root = scipy.optimize.brentq(
        bessel_function_of_nu, bracket_interval[0], bracket_interval[1]
    )

    # The relationship between x and the found root nu is: x = nu + 1.
    largest_x_value = largest_nu_root + 1

    # Output the steps of the solution
    # The 'final equation' is I_nu(2)=0, where nu and x are the numbers.
    print(f"The summation is equal to the modified Bessel function I_nu(2) with nu = x - 1.")
    print(f"We solve the equation I_nu(2) = 0 for the largest root nu.")
    print(f"Largest root nu found: {largest_nu_root}")
    print(f"Corresponding x value (x = nu + 1): {largest_x_value}")

    # Print the final answer in the required format {-a.bbb}
    print("\nFinal Answer:")
    print(f"{{{largest_x_value:.3f}}}")

except (ImportError, ValueError) as e:
    print(f"An error occurred: {e}")
    print("Please ensure you have the scipy library installed (`pip install scipy`).")
    print("This script cannot be run without it.")
