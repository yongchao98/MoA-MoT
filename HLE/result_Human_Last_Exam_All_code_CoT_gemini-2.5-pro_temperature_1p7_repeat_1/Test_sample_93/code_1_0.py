import math
from collections import Counter

def solve_sum_of_squares():
    """
    Finds the number of non-negative integer solutions to x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = 2024.
    """
    target = 2024
    num_vars = 5
    limit = int(math.sqrt(target))
    
    unique_solutions = []
    
    for x1 in range(limit + 1):
        s1 = x1 * x1
        if s1 > target:
            break
        for x2 in range(x1 + 1):
            s2 = s1 + x2 * x2
            if s2 > target:
                break
            for x3 in range(x2 + 1):
                s3 = s2 + x3 * x3
                if s3 > target:
                    break
                for x4 in range(x3 + 1):
                    s4 = s3 + x4 * x4
                    if s4 > target:
                        break
                    
                    rem = target - s4
                    if rem < 0:
                        continue
                        
                    x5 = int(math.sqrt(rem))
                    
                    if x5 * x5 == rem:
                        if x5 <= x4:
                            solution = (x1, x2, x3, x4, x5)
                            unique_solutions.append(solution)

    total_solutions = 0
    print(f"Finding non-negative integer solutions for x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = {target}:\n")
    
    if not unique_solutions:
        print("No solutions found.")
    else:
        print("Unique solution sets (in non-increasing order) and their contributions:")
        for s in unique_solutions:
            counts = Counter(s)
            
            numerator = math.factorial(num_vars)
            denominator = 1
            for count in counts.values():
                denominator *= math.factorial(count)
            
            permutations = numerator // denominator
            total_solutions += permutations
            
            # Print the equation for each solution
            equation_str_parts = [f"{val}^2" for val in s]
            equation_str = " + ".join(equation_str_parts)
            
            values_str_parts = [f"{val*val}" for val in s]
            values_str = " + ".join(values_str_parts)

            print(f"Solution: ({', '.join(map(str, s))})")
            print(f"  {equation_str} = {values_str} = {target}")
            print(f"  Number of permutations: {permutations}\n")
    
    print("-" * 30)
    print(f"Total number of non-negative integer solutions: {total_solutions}")

solve_sum_of_squares()
<<<4560>>>