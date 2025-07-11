import math

def solve_random_walk_probability():
    """
    Calculates the probability for a 2D random walk to visit a target set
    before leaving a disk, using a continuous approximation (Brownian motion).
    """
    # Parameters from the problem
    R = 1000.0  # Radius of the outer disk
    
    # Starting point S = (0, 300)
    start_point = (0, 300)
    
    # Target set A = {(0,0), (2,0)}
    # The center of the target set is (1,0)
    target_center = (1, 0)
    
    # The distance 'd' from the starting point to the center of the target set
    d = math.sqrt((start_point[0] - target_center[0])**2 + (start_point[1] - target_center[1])**2)
    
    # The effective radius 'epsilon' of the target set.
    # Approximating the two points by a line segment of length L=2.
    # The logarithmic capacity is L/4, which we use as the effective radius.
    L = 2.0
    epsilon = L / 4.0
    
    # The probability P is given by the formula for hitting a small disk before a large one.
    # P = ln(R / d) / ln(R / epsilon)
    # The base of the logarithm does not matter as long as it is consistent.
    probability = math.log(R / d) / math.log(R / epsilon)
    
    # Print the equation and the values used
    print("The probability P is calculated using the formula for hitting probability in Brownian motion:")
    print("P = ln(R / d) / ln(R / epsilon)")
    print(f"where:")
    print(f"R (Radius of boundary disk) = {R}")
    print(f"d (Distance from start to target center) = sqrt({(start_point[0] - target_center[0])**2} + {start_point[1] - target_center[1]}**2) = {d:.4f}")
    print(f"epsilon (Effective radius of target set) = {epsilon}")
    print("\nPlugging in the values:")
    print(f"P = ln({R} / {d:.4f}) / ln({R} / {epsilon})")
    print(f"P = ln({R/d:.4f}) / ln({R/epsilon:.4f})")
    print(f"P = {math.log(R / d):.4f} / {math.log(R / epsilon):.4f}")
    print(f"P = {probability:.3f}")

solve_random_walk_probability()