import math
from collections import Counter

def solve_diophantine_equation():
    """
    Finds the number of non-negative integer solutions to the equation:
    x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024
    """
    target = 2024
    limit = int(math.sqrt(target))
    
    # Pre-calculate squares for performance improvement
    squares = [i*i for i in range(limit + 1)]
    
    # Use a set to store unique unordered solutions (as sorted tuples)
    solution_sets = set()

    # Find solutions with the constraint x1 >= x2 >= x3 >= x4 >= x5 >= 0
    # This avoids finding duplicate sets and speeds up the search.
    for i1 in range(limit, -1, -1):
        s1 = squares[i1]
        
        limit2 = int(math.sqrt(target - s1))
        for i2 in range(min(i1, limit2), -1, -1):
            s2 = squares[i2]
            
            limit3 = int(math.sqrt(target - s1 - s2))
            for i3 in range(min(i2, limit3), -1, -1):
                s3 = squares[i3]

                limit4 = int(math.sqrt(target - s1 - s2 - s3))
                for i4 in range(min(i3, limit4), -1, -1):
                    s4 = squares[i4]
                    
                    remaining_sum = target - (s1 + s2 + s3 + s4)
                    
                    # The value for x5^2 must be non-negative.
                    if remaining_sum >= 0:
                        i5 = int(math.sqrt(remaining_sum))
                        
                        # Check if remaining_sum is a perfect square and ordering is maintained
                        if squares[i5] == remaining_sum and i5 <= i4:
                            solution_tuple = (i1, i2, i3, i4, i5)
                            solution_sets.add(solution_tuple)

    total_solutions = 0
    fact5 = math.factorial(5)

    print(f"Finding non-negative integer solutions for x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = {target}\n")
    print("Found the following unique solution sets (and their permutations):\n")

    # Sort sets for deterministic output
    sorted_sets = sorted(list(solution_sets), reverse=True)

    for sol_set in sorted_sets:
        # Count frequency of each number in the set
        counts = Counter(sol_set)
        denominator = 1
        for count in counts.values():
            denominator *= math.factorial(count)
        
        num_perms = fact5 // denominator
        total_solutions += num_perms
        
        # Output the numbers in the equation
        squares_str = [f"{x}^2" for x in sol_set]
        values_str = [f"{x*x}" for x in sol_set]
        
        print(f"Solution Set: {sol_set}")
        print(f"Equation: {' + '.join(squares_str)} = {' + '.join(values_str)} = {target}")
        print(f"Number of distinct permutations: {num_perms}\n")

    print("----------------------------------------------------")
    print(f"Total number of non-negative integer solutions: {total_solutions}")
    return total_solutions

if __name__ == '__main__':
    final_answer = solve_diophantine_equation()
    print(f"\n<<< {final_answer} >>>")
