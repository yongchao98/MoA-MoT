import math
from collections import Counter

def solve_diophantine_equation():
    """
    This function finds the number of non-negative integer solutions to the equation:
    x₁² + x₂² + x₃² + x₄² + x₅² = 2024
    """
    N = 2024
    limit = math.isqrt(N)
    total_solutions = 0
    base_solutions = []

    # Iterate to find unique ordered sets (multisets) of solutions {x₁, ..., x₅}
    # where x₁ ≤ x₂ ≤ x₃ ≤ x₄ ≤ x₅. This ensures each set is found only once.
    for x1 in range(limit + 1):
        s1 = x1 * x1
        for x2 in range(x1, limit + 1):
            s2 = x2 * x2
            sum2 = s1 + s2
            if sum2 > N: break
            for x3 in range(x2, limit + 1):
                s3 = x3 * x3
                sum3 = sum2 + s3
                if sum3 > N: break
                for x4 in range(x3, limit + 1):
                    s4 = x4 * x4
                    sum4 = sum3 + s4
                    if sum4 > N: break
                    
                    remaining = N - sum4
                    # Efficiently check if 'remaining' is a perfect square
                    x5 = math.isqrt(remaining)
                    
                    if x5 * x5 == remaining and x5 >= x4:
                        solution_set = (x1, x2, x3, x4, x5)
                        base_solutions.append(solution_set)

    print(f"Finding non-negative integer solutions for the equation:\n"
          f"x₁² + x₂² + x₃² + x₄² + x₅² = {N}\n")

    # Calculate total solutions by summing permutations of each unique set
    for sol in base_solutions:
        counts = Counter(sol)
        # Calculate permutations for the multiset: 5! / (count1! * count2! * ...)
        perms = math.factorial(5)
        for count in counts.values():
            perms //= math.factorial(count)
        
        total_solutions += perms
        
        # Output the equation for each base solution
        equation_parts = [f"{x}²" for x in sol]
        sum_parts = [f"{x*x}" for x in sol]
        print(f"Base Solution: {sol}")
        print(f"  {' + '.join(equation_parts)} = {' + '.join(sum_parts)} = {N}")
        print(f"  Number of permutations for this set: {perms}\n")
    
    print("-" * 40)
    print(f"Found {len(base_solutions)} unique sets of integers (base solutions).")
    print(f"The total number of non-negative integer solutions is {total_solutions}.")

# Execute the function to find and print the solution
solve_diophantine_equation()