import numpy as np
from scipy.integrate import quad

def borwein_integrand(x, n):
    """
    The integrand of the Borwein integral I_n.
    The function is a product of n sinc functions.
    sinc(x) = sin(x)/x
    """
    # Handle the limit at x -> 0 where sinc(x) is 1.
    if x == 0:
        return 1.0
    
    product = 1.0
    for k in range(1, n + 1):
        arg = x / k
        product *= np.sin(arg) / arg
    return product

def analyze_borwein_integrals(max_n=9):
    """
    Calculates and analyzes Borwein integrals I_n for n=1 to max_n.
    """
    print("--- Analysis of Borwein Integrals I_n ---")
    
    # Value of pi/2 for comparison
    pi_half = np.pi / 2
    
    # To store previous value for monotonicity check
    last_I_n = None 

    for n in range(1, max_n + 1):
        print(f"\n--- n = {n} ---")
        
        # Perform numerical integration from 0 to infinity
        # The quad function returns the integral value and an estimated error.
        val, err = quad(borwein_integrand, 0, np.inf, args=(n,))
        
        print(f"Calculated I_{n} = {val:.14f}")
        print(f"pi/2             = {pi_half:.14f}")
        
        # Check proposition P(n): I_n = pi/2
        # We check if the difference is within a very small tolerance.
        is_p_n_true = np.isclose(val, pi_half)
        diff = val - pi_half
        
        print(f"Proposition P({n}) [I_{n} = pi/2] is: {is_p_n_true} (Difference: {diff:.2e})")

        # Check statement C: If P(n) is false, then I_n < pi/2
        if not is_p_n_true:
            print(f"Check for statement C: I_{n} < pi/2 is {val < pi_half}")

        # Check statement F for n=5
        if n == 5:
            is_f_true = abs(diff) < 1e-5
            print(f"Check for statement F: |I_5 - pi/2| < 10^-5 is {is_f_true}")

        # Check statement G: The sequence {I_n} is monotonically decreasing (I_n <= I_{n-1})
        if last_I_n is not None:
            # We use isclose for the equality part
            is_monotonic = (val < last_I_n) or np.isclose(val, last_I_n)
            print(f"Check for statement G: I_{n} <= I_{n-1} is {is_monotonic}")

        last_I_n = val

analyze_borwein_integrals()