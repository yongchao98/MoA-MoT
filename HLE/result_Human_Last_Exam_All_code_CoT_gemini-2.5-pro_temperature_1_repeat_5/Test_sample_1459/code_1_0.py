import math

def solve_gh_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle with the intrinsic metric by solving a transcendental equation
    numerically using the bisection method.
    """
    # The Gromov-Hausdorff distance 'd' is the solution to the equation:
    # 1 = 2 * d * arccos(2/d - 1)
    # We are looking for the root of the function g(d) = 2 * d * arccos(2/d - 1) - 1
    
    def g(d):
        # The argument of arccos must be in [-1, 1].
        # 2/d - 1 is in [-1, 1] if d >= 1.
        if d < 1:
            return float('nan') # Not in the domain of the solution
        arg = 2 / d - 1
        # Handle potential floating point inaccuracies at the boundary
        if arg > 1.0:
            arg = 1.0
        if arg < -1.0:
            arg = -1.0
        return 2 * d * math.acos(arg) - 1

    # We know the solution lies in the interval [1, 2] because:
    # g(1) = -1
    # g(2) = 2*pi - 1 > 0
    a, b = 1.0, 2.0

    # Bisection method parameters
    tolerance = 1e-12
    max_iterations = 100

    g_a = g(a)

    for _ in range(max_iterations):
        mid = (a + b) / 2
        if b - a < tolerance:
            break
        
        g_mid = g(mid)
        
        if g_mid == 0:
            a = b = mid
            break
        
        if math.copysign(1, g_mid) == math.copysign(1, g_a):
            a = mid
            g_a = g_mid
        else:
            b = mid
            
    distance = (a + b) / 2
    
    # Output the explanation and the result
    print("The problem is to find the Gromov-Hausdorff distance, d, between the interval [0,1] and the unit circle.")
    print("This requires solving the following transcendental equation for d:")
    
    # Outputting each number in the final equation as requested.
    L = 1
    D = 2
    print(f"{L} = 2 * d * arccos({D} / d - 1)")
    
    print("\nSolving this equation numerically gives the distance:")
    print(f"{distance:.10f}")


if __name__ == '__main__':
    solve_gh_distance()
