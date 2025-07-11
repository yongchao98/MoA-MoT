import math

def solve():
    """
    Calculates the probability that a 2D random walk, conditioned to never enter the origin,
    starting from (0,1), will eventually hit the set of neighbors of (3600,0).

    The probability P is given by the ratio of potential kernel values: a(x0) / a(y),
    where x0 = (0,1) is the starting point and y = (3600,0) is the target point.

    The potential kernel a(x) for large |x| is approximated by:
    a(x) ~= (2/pi) * ln|x| + K
    where K = (2/pi) * (gamma + 1.5 * ln(2)), and gamma is the Euler-Mascheroni constant.
    For x0=(0,1), |x0|=1, so a(x0) is approximated by K.
    """

    # Constants
    gamma = 0.5772156649
    pi = math.pi
    ln2 = math.log(2)

    # Coordinates
    x0_norm = 1.0
    y_norm = 3600.0

    # Calculate the constant K in the asymptotic expansion of the potential kernel
    # K is a good approximation for a(1,0)
    K = (2 / pi) * (gamma + 1.5 * ln2)
    a_x0 = K

    # Calculate the potential kernel value at y = (3600,0)
    ln_y_norm = math.log(y_norm)
    a_y = (2 / pi) * ln_y_norm + K
    
    # The probability is the ratio of the potential kernel values
    probability = a_x0 / a_y
    
    # Print the equation with calculated values
    print("The probability P is approximated by the ratio a(x0) / a(y).")
    print("x0 = (0,1), y = (3600,0)")
    print(f"a(x) is approximated by (2/pi)*ln|x| + K for large |x|.")
    print(f"K = (2/pi)*(gamma + 1.5*ln(2))")
    print(f"  ~= (2/{pi:.4f})*({gamma:.4f} + 1.5*{ln2:.4f}) ~= {K:.4f}")
    print(f"a(x0) = a(0,1) is approximated by K ~= {a_x0:.4f}")
    print(f"a(y) = a(3600,0) is approximated by (2/pi)*ln(3600) + K")
    print(f"  ~= (2/{pi:.4f})*{ln_y_norm:.4f} + {K:.4f} ~= {a_y:.4f}")
    print("\nFinal calculation:")
    print(f"P ~= a(x0) / a(y) = {a_x0:.4f} / {a_y:.4f} = {probability:.4f}")
    
    # Output the result with two significant digits
    print(f"\nThe approximate answer with two significant digits is: {probability:.2g}")

solve()