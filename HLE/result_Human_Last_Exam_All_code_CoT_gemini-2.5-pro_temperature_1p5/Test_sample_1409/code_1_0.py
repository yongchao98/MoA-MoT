import math

def solve_prism_permutations():
    """
    Finds the number of permutations of dimensions (x, y, z) for a rectangular prism
    with volume 216, satisfying the condition x^(log_2(y)) = y^(log_2(z)).
    """
    volume = 216
    
    # The logarithmic condition simplifies to (y = 1) or (x = z).
    # These two cases are disjoint because if y=1 and x=z, then x^2 = 216, which has no integer solution for x.
    # So we can count the solutions for each case and add them up.

    print("The problem requires finding the number of ordered triples (x, y, z) of positive integers where:")
    print("1. x * y * z = 216")
    print("2. x^(log₂(y)) = y^(log₂(z)) which simplifies to (y = 1) or (x = z)")
    print("-" * 50)

    # Case 1: y = 1
    # The volume equation becomes x * z = 216.
    # We need to find all pairs of factors (x, z) for 216.
    case1_solutions = []
    for x in range(1, volume + 1):
        if volume % x == 0:
            z = volume // x
            case1_solutions.append((x, 1, z))
            
    print("Solutions for Case 1 (y = 1):")
    for sol in case1_solutions:
        print(sol)
    print(f"\nNumber of permutations in Case 1: {len(case1_solutions)}")
    print("-" * 50)
    
    # Case 2: x = z
    # The volume equation becomes x^2 * y = 216.
    # We need to find all x such that x^2 is a divisor of 216.
    case2_solutions = []
    # We only need to check x up to the square root of the volume.
    limit = int(math.sqrt(volume))
    for x in range(1, limit + 1):
        if (x * x > 0) and (volume % (x * x) == 0):
            y = volume // (x * x)
            case2_solutions.append((x, y, x))
            
    print("Solutions for Case 2 (x = z):")
    for sol in case2_solutions:
        print(sol)
    print(f"\nNumber of permutations in Case 2: {len(case2_solutions)}")
    print("-" * 50)
    
    # The total number of permutations is the sum of the counts from the two disjoint cases.
    total_permutations = len(case1_solutions) + len(case2_solutions)
    
    print("The total number of permutations is the sum of the counts from the two disjoint cases.")
    print(f"Total = {len(case1_solutions)} (from y=1) + {len(case2_solutions)} (from x=z)")
    print(f"Total = {total_permutations}")

solve_prism_permutations()