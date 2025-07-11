import math

def find_t4():
    """
    Finds T(4), the sum of all elements in S(4).
    S(n) is the set of all numbers that can be expressed as a sum of n distinct
    positive integers whose reciprocals sum to exactly 1.
    """
    solution_sets = []
    
    # Based on the analysis 1 < x1 < 4, so x1 can be 2 or 3.
    # We iterate through possible x1 and x2 to find all solution sets.
    for x1 in range(2, 4):
        # We derive the upper bound for x2 from 1 - 1/x1 < 3/x2
        x2_upper_bound = math.ceil(3 * x1 / (x1 - 1))
        for x2 in range(x1 + 1, x2_upper_bound):
            
            # The equation for x3 and x4 is 1/x3 + 1/x4 = 1 - 1/x1 - 1/x2
            # Let's express the right side as a fraction A/B
            # A = x1*x2 - x1 - x2
            # B = x1*x2
            # This can be transformed into (A*x3 - B) * (A*x4 - B) = B*B
            
            A = x1 * x2 - x1 - x2
            if A <= 0:  # This happens if 1/x1 + 1/x2 >= 1, no positive solutions for x3, x4
                continue
            
            B = x1 * x2
            B_squared = B * B

            # We look for factor pairs (u, v) of B_squared to solve for x3 and x4
            for u in range(1, int(math.sqrt(B_squared)) + 1):
                if B_squared % u == 0:
                    v = B_squared // u
                    
                    # From u = A*x3 - B, we get x3 = (B+u)/A
                    # From v = A*x4 - B, we get x4 = (B+v)/A
                    if (B + u) % A == 0 and (B + v) % A == 0:
                        x3 = (B + u) // A
                        x4 = (B + v) // A
                        
                        # We require distinct integers x1 < x2 < x3 < x4
                        if x3 > x2 and x3 < x4:
                            solution_sets.append([x1, x2, x3, x4])

    # S(4) is the set of the sums of the numbers in each solution set.
    s4 = set()
    for s in solution_sets:
        s4.add(sum(s))
        
    # T(4) is the sum of all elements in S(4).
    t4 = sum(s4)
    
    sorted_sums = sorted(list(s4))
    
    # Print the final equation as requested.
    sum_str = " + ".join(map(str, sorted_sums))
    print(f"The elements of S(4) are the sums of each valid set of integers.")
    print(f"The sets are: {solution_sets}")
    print(f"The sums (S(4)) are: {sorted_sums}")
    print("\nThe final equation for T(4) is:")
    print(f"{sum_str} = {t4}")

find_t4()
<<<208>>>