import math

def solve_polygon_area_ratio(n):
    """
    Calculates the area ratio of an n-sided polygon constructed by extending
    the alternate edges of a 2n-sided regular polygon.
    """
    if n < 3:
        print("A regular polygon must have at least 3 sides.")
        return
        
    print(f"This script calculates the area ratio for n = {n}.")
    print("The general formula for the ratio is: (1 + sec(pi/n)) / 2")
    print("-" * 50)
    
    # --- Step 1: Calculate pi / n ---
    pi_over_n = math.pi / n
    
    # --- Step 2: Calculate sec(pi/n) ---
    # We calculate cos(pi/n) first, as sec(x) = 1/cos(x)
    cos_pi_over_n = math.cos(pi_over_n)
    if cos_pi_over_n == 0:
        print("Calculation Error: cos(pi/n) is zero, secant is undefined.")
        return
    sec_pi_over_n = 1 / cos_pi_over_n
    
    # --- Step 3: Calculate the final ratio using the derived values ---
    # The final equation is (1 + sec_pi_over_n) / 2
    numerator = 1 + sec_pi_over_n
    denominator = 2
    ratio = numerator / denominator

    print("The numbers that go into the final equation are:")
    print(f"1 (a constant)")
    print(f"sec(pi/{n}) = {sec_pi_over_n}")
    print(f"2 (a constant)")
    print("")

    print(f"The final equation with the calculated values is:")
    print(f"Ratio = ({1} + {sec_pi_over_n}) / {denominator}")
    print(f"\nThe area of the {n}-sided polygon is {ratio} times larger than the {2*n}-sided polygon.")


# Using n=3 as given in the problem statement (triangle from a hexagon).
# The expected result is 3/2 = 1.5.
solve_polygon_area_ratio(3)
