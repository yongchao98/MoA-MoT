import numpy as np

def F(X, Y):
    """
    This function represents the equation of the curve Gamma in the transformed coordinates X=x+y and Y=x-y.
    A point is inside the curve if F(X,Y) < 0.
    """
    X = float(X)
    Y = float(Y)
    term1 = 1200.0 * Y**4
    term2 = -20.0 * X**3 * Y**2
    term3 = 3.0 * (X**2 - 400.0)**3
    return term1 + term2 + term3

def solve():
    """
    Counts the number of poles of f(z) inside the contour Gamma.
    """
    pole_count = 0
    a_min = -2024
    a_max = 2024

    # Based on manual analysis, poles only exist for |n| <= 3.
    # We check a slightly larger range to be safe.
    n_min = -5
    n_max = 5

    for n in range(n_min, n_max + 1):
        for a in range(a_min, a_max + 1):
            # Coordinates of the pole
            x = a
            y = 2 * np.pi * n
            
            # Transformed coordinates
            X = x + y
            Y = x - y
            
            # Check if the pole is inside the curve
            if F(X, Y) < 0:
                pole_count += 1
    
    # The integral is pole_count * 2 * pi * i
    # The question requires printing each number in the final equation.
    # The final equation is Integral = N * 2 * pi * i, where N is the number of poles.
    print(f"The number of poles inside the contour is: {pole_count}")
    print("The value of the contour integral is given by the equation:")
    print(f"oint_Gamma f(z) dz = {pole_count} * 2 * pi * i")

if __name__ == '__main__':
    solve()