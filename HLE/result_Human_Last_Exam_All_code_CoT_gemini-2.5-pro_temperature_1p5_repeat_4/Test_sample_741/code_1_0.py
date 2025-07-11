import scipy.special as sp
import scipy.optimize as opt
import numpy as np

# This script requires the scipy library. You can install it using:
# pip install scipy

# The problem is to find the largest x for which the sum is zero.
# The sum is equivalent to the modified Bessel function of the first kind, I_v(z),
# where v = x - 1 and z = 2.
# So we need to solve I_{x-1}(2) = 0.
# Let v = x - 1. We need to find the largest root of I_v(2) = 0.

# Define the function whose root we want to find.
def bessel_function_of_order_v(v):
  """
  Represents the modified Bessel function I_v(2).
  """
  return sp.iv(v, 2.0)

try:
  # The roots of I_v(z) for v are all negative.
  # The largest root (closest to zero) is known to be in the interval [-2, -1].
  # We use the Brent-Dekker method to find the root in this bracket.
  v_root = opt.brentq(bessel_function_of_order_v, -2.0, -1.0)
  
  # Now, we calculate x from v using x = v + 1.
  x_root = v_root + 1
  
  # Print the numbers used in the final calculation, as requested.
  print(f"The largest root for v in I_v(2)=0 is found to be: {v_root}")
  print(f"The value for x is calculated from v using x = v + 1.")
  print(f"Final equation: {v_root:.3f} + 1 = {x_root:.3f}")

  # Print the final answer in the specified format {-a.bbb}
  print("\nFormatted answer:")
  # Use numpy.format_float_scientific to handle the negative zero case gracefully
  # and ensure the format {-0.xyz} is printed correctly.
  formatted_x = np.format_float_positional(x_root, precision=3, fractional=False, trim='-')
  print(f"{{-{formatted_x}}}")


except ImportError:
  print("SciPy library not found. Please install it using 'pip install scipy'")
except ValueError as e:
  print(f"Root finding failed. The function may not have a root in the given interval [-2.0, -1.0]. Error: {e}")
except Exception as e:
  print(f"An unexpected error occurred: {e}")
