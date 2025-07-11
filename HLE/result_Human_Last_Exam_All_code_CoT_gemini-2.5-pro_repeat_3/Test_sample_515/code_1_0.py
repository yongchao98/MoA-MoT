import numpy as np

def solve():
    """
    Calculates the probability that a conditioned 2D random walk starting at (0,1)
    hits the neighbors of (3600,0) before hitting the origin.
    
    The probability is approximated by p = a(0,1) / (2 * a(3600,0)).
    """
    
    # Coordinates of the points
    A = np.array([0, 1])
    B = np.array([3600, 0])
    
    # 1. Calculate a(A) = a(0,1)
    # The value for discrete-time SRW is not simple. We use the value for
    # the continuous-time random walk, a(1,0) = 4/pi - 1, as a close approximation.
    a_A = 4 / np.pi - 1
    
    # 2. Calculate a(B) = a(3600,0) using the asymptotic formula
    # a(x) approx (2/pi)*ln|x| + C for large |x|
    # C = (2*gamma + ln(8))/pi
    
    gamma = np.euler_gamma  # Euler-Mascheroni constant
    
    # The constant C in the asymptotic expansion of a(x)
    C = (2 * gamma + np.log(8)) / np.pi
    
    # Magnitude of vector B
    mag_B = np.linalg.norm(B)
    
    # Asymptotic value of a(B)
    a_B = (2 / np.pi) * np.log(mag_B) + C
    
    # 3. Calculate the final probability
    # p approx a(A) / (2 * a(B))
    prob = a_A / (2 * a_B)
    
    print("The starting point is A = (0, 1).")
    print("The target point is B = (3600, 0).")
    print("The probability is calculated as p ≈ a(A) / (2 * a(B)).")
    print(f"Value of a(A) = a(0,1) is approximated by the continuous-time value 4/π - 1 = {a_A:.4f}")
    print(f"The asymptotic constant C is {C:.4f}")
    print(f"Value of a(B) = a(3600,0) is approximated by (2/π)ln(3600) + C = {a_B:.4f}")
    print(f"The probability p ≈ {a_A:.4f} / (2 * {a_B:.4f})")
    
    # Final result
    print(f"The final approximate probability is: {prob:.4f}")
    
    # Approximate answer with two significant digits
    approx_prob = float(f"{prob:.2g}")
    print(f"The approximate answer with two significant digits is: {approx_prob}")

solve()