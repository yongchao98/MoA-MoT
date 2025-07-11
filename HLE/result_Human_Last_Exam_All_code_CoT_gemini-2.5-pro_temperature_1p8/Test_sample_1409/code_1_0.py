import math

def solve_prism_permutations():
    """
    Finds the number of permutations of dimensions (x, y, z) for a rectangular prism
    with volume 216 and satisfying the condition x^(log_2 y) = y^(log_2 z).
    """
    volume = 216
    solutions = set()

    # Case 1: y = 1
    # The condition simplifies to x * z = 216.
    # We find all divisors of 216 for x, which determines z.
    case1_solutions = []
    for x in range(1, volume + 1):
        if volume % x == 0:
            z = volume // x
            # y is fixed at 1
            y = 1
            solution = (x, y, z)
            case1_solutions.append(solution)
            solutions.add(solution)
    
    count_case1 = len(case1_solutions)

    # Case 2: x = z
    # The condition simplifies to x^2 * y = 216.
    # We find all x where x^2 is a divisor of 216.
    case2_solutions = []
    # We only need to check x up to the square root of the volume.
    for x in range(1, int(math.sqrt(volume)) + 1):
        if volume % (x*x) == 0:
            y = volume // (x*x)
            # z is equal to x
            z = x
            solution = (x, y, z)
            case2_solutions.append(solution)
            solutions.add(solution)

    count_case2 = len(case2_solutions)

    print(f"The problem is solved by finding all integer triples (x, y, z) where x*y*z = 216 and either y=1 or x=z.")
    print("-" * 50)
    
    # "output each number in the final equation"
    print(f"Number of solutions where y=1: {count_case1}")
    print(f"Number of solutions where x=z: {count_case2}")
    
    # As reasoned in the thinking steps, the cases are disjoint, so a simple sum is sufficient.
    # The use of a set `solutions` also correctly calculates the union of the two sets of results.
    total_permutations = len(solutions)
    
    print("-" * 50)
    print(f"The final equation for the total number of permutations is the sum of the counts from the two disjoint cases:")
    print(f"{count_case1} + {count_case2} = {total_permutations}")
    
solve_prism_permutations()
