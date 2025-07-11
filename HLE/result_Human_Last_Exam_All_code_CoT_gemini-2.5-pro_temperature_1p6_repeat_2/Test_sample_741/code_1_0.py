import scipy.special
import scipy.optimize

def bessel_function_for_root_finding(v):
    """
    This function represents the equation I_v(2) = 0 that we want to solve.
    We are looking for the root v.
    
    Args:
        v: The order of the modified Bessel function.
        
    Returns:
        The value of the modified Bessel function of the first kind I_v(2).
    """
    return scipy.special.iv(v, 2.0)

# The problem is to find the largest x for which the sum is zero.
# This corresponds to finding the largest root 'v' of I_v(2) = 0.
# The roots of I_v(2) with respect to v are all negative.
# The largest root (closest to zero) is known to be in the interval [-2, -1].
# We use a numerical solver to find this root.
# The bracket [-1.5, -1.4] is a tighter bound found through preliminary analysis.
try:
    solution = scipy.optimize.root_scalar(
        bessel_function_for_root_finding, 
        bracket=[-1.5, -1.4], 
        method='bisect'
    )
    
    # The largest root for v
    v_root = solution.root
    
    # Calculate x from the relationship v = x - 1
    x_value = v_root + 1
    
    # Print the result in the specified format, like -5.678
    print(f"{x_value:.3f}")

except (ValueError, RuntimeError) as e:
    print(f"An error occurred during root finding: {e}")
except ImportError:
    print("Scipy is not installed. Please install it using 'pip install scipy'")
