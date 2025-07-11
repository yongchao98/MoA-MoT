import math
from collections import Counter
import itertools

def solve_diophantine_equation():
    """
    This function finds the number of non-negative integer solutions to the equation:
    x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024
    """

    target = 2024
    # The maximum possible value for any x_i is floor(sqrt(2024)) = 44
    limit = int(math.sqrt(target))
    
    ordered_solution_sets = []

    # We generate combinations of 5 numbers from 0 to 44 with replacement.
    # This ensures that we find all unique sets of numbers {x1, x2, x3, x4, x5}
    # where the order doesn't matter yet, and repetitions are allowed.
    # `itertools.combinations_with_replacement` generates them in sorted (non-decreasing) order.
    possible_values = range(limit + 1)
    
    for combo in itertools.combinations_with_replacement(possible_values, 5):
        current_sum_sq = sum(x*x for x in combo)
        if current_sum_sq == target:
            # We found a set of numbers that works.
            # Store it as a tuple in our list of solutions.
            ordered_solution_sets.append(combo)

    # Now, for each unique set of solutions, we calculate the number of distinct
    # permutations. For a set {x1, x2, x3, x4, x5}, the number of permutations is
    # 5! / (n1! * n2! * ...), where n_k is the frequency of each distinct number in the set.
    
    total_solutions = 0
    print(f"Found {len(ordered_solution_sets)} unique sets of non-negative integers (unordered) whose squares sum to {target}.")
    print("----------------------------------------------------------------------")
    
    if not ordered_solution_sets:
        print("No solutions found.")
    else:
        # Using math.factorial is cleaner
        for solution_set in ordered_solution_sets:
            # Calculate number of permutations for this set
            counts = Counter(solution_set)
            denominator = 1
            for count in counts.values():
                denominator *= math.factorial(count)
            
            num_permutations = math.factorial(5) // denominator
            total_solutions += num_permutations

            # Format the equation string as requested
            equation_str_parts = [f"{x}^2" for x in solution_set]
            equation_sum_parts = [f"{x*x}" for x in solution_set]

            print(f"Solution set: {solution_set}")
            print(f"Equation: {' + '.join(equation_str_parts)} = {' + '.join(equation_sum_parts)} = {target}")
            print(f"Number of distinct permutations (solutions): {num_permutations}")
            print("----------------------------------------------------------------------")

    print(f"\nThe total number of non-negative integer solutions is: {total_solutions}")

solve_diophantine_equation()
<<<930>>>