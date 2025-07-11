import math

def calculate_values(n):
    """
    This function performs the calculation that determines the solution.
    
    After a conformal map z -> 1/z:
    - The original target center (n, 0) maps to approximately (1/n, 0).
    - The original target radius n^(1/3) maps to a new radius of n^(-5/3).
    
    The nearest integer lattice point to the new center is the origin (0,0).
    Its distance from the new center is 1/n.
    
    The random walk can only hit the target if the target disk contains an integer lattice point.
    This requires the distance to the nearest point to be less than the radius.
    """
    if n <= 1:
        print(f"The analysis holds for large n (n>1).")
        return

    # Squared distance from the mapped center (1/n, 0) to the nearest lattice point (0,0)
    dist_sq = (1/n)**2
    
    # Squared radius of the mapped target disk
    radius_sq = (n**(-5/3))**2
    
    print(f"For n = {n}:")
    print(f"The nearest lattice point is (0,0), at distance 1/n from the center of the mapped target.")
    print(f"Equation: is (distance)^2 < (radius)^2 ?")
    # We output each number in the equation, as requested.
    print(f"Is {dist_sq:.2e} < {radius_sq:.2e} ?")
    
    # Check the inequality
    contains_point = dist_sq < radius_sq
    print(f"Result: {contains_point}\n")
    
    # Analytical conclusion for n > 1:
    # We are comparing 1/n with n**(-5/3).
    # This is equivalent to comparing n**(-1) with n**(-5/3).
    # Since the exponent -1 is greater than -5/3, for n>1, n**(-1) > n**(-5/3).
    # The distance is always greater than the radius for n > 1.
    

# Demonstrate the calculation for a large n
calculate_values(100)

# The analysis shows that for any n > 1, the mapped target disk contains no
# lattice points. Therefore, the probability p_n is 0 for n >= 2.
# The limit must therefore be 0.
final_limit = 0

print("Final Conclusion:")
print("lim_{n->inf} p_n = 0")
