import math

def solve_prism_permutations():
    """
    This function finds the number of ordered triplets (x, y, z) that satisfy
    the given conditions for a rectangular prism.
    """
    volume = 216
    
    print("The problem is to find the number of ordered integer triplets (x, y, z) such that:")
    print(f"1. x * y * z = {volume}")
    print("2. x^(log₂(y)) = y^(log₂(z))")
    print("\nThe second condition simplifies to (y=1) or (x=z).")
    print("We count the solutions for each case. The total is the sum of these counts since their intersection is empty for integers.\n")

    # Case 1: y = 1
    # This implies x * z = 216. We count the number of divisors of 216.
    count_case1 = 0
    for i in range(1, volume + 1):
        if volume % i == 0:
            count_case1 += 1
            
    print(f"Case 1: Find solutions where y = 1.")
    print(f"This requires x * z = {volume}. The number of such integer pairs (x, z) is the number of divisors of {volume}.")
    print(f"The number of solutions for this case is {count_case1}.")

    # Case 2: x = z
    # This implies x^2 * y = 216.
    count_case2 = 0
    # We check for x values where x*x divides the volume.
    for x in range(1, int(math.sqrt(volume)) + 1):
        if (x * x) != 0 and volume % (x * x) == 0:
            count_case2 += 1
            
    print(f"\nCase 2: Find solutions where x = z.")
    print(f"This requires x² * y = {volume}. This holds true for every x where x² is a perfect square divisor of {volume}.")
    print(f"The number of solutions for this case is {count_case2}.")
    
    # Total number of permutations
    total_solutions = count_case1 + count_case2
    
    print("\nSince there are no integer solutions where y=1 and x=z simultaneously, the two sets of solutions are disjoint.")
    print("The final calculation is:")
    print(f"{count_case1} + {count_case2} = {total_solutions}")

solve_prism_permutations()