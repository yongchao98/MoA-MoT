import scipy.special
import scipy.optimize

# The given summation can be expressed in terms of the modified Bessel function
# of the first kind, I_v(z), as S(x) = I_{x-1}(2).
# We need to find the largest x for which S(x) = 0.
# This is equivalent to finding the largest root of I_v(2) = 0, where v = x - 1.

def bessel_function_of_order_v(v):
  """
  Represents the function I_v(2) for which we want to find the root.
  The input 'v' is the order of the Bessel function.
  """
  return scipy.special.iv(v, 2)

# The real roots 'v' for I_v(2)=0 are all negative. We need to find the largest one.
# By evaluating the function at different points, we can locate an interval for the largest root.
# For example, iv(-2.2, 2) is positive, and iv(-2.3, 2) is negative.
# So, the largest root lies in the interval [-2.3, -2.2].
lower_bound_v = -2.3
upper_bound_v = -2.2

# Use the brentq method to find the root v in the specified interval.
try:
    v_root = scipy.optimize.brentq(bessel_function_of_order_v, lower_bound_v, upper_bound_v)

    # The value of x is derived from the relation x = v + 1.
    x_root = v_root + 1

    # The final equation is x = v + 1. We print each number involved.
    print(f"The largest root for the order v is: {v_root:.10f}")
    print(f"The constant in the equation is: 1")
    print(f"The largest value for x is: {x_root:.10f}")

    # Print the result in the requested format {-a.bbb}
    print(f"\nFormatted Answer: {{{x_root:.3f}}}")

except (ValueError, RuntimeError) as e:
    print(f"An error occurred during root finding: {e}")
except ImportError:
    print("SciPy library is required. Please install it using 'pip install scipy'")
