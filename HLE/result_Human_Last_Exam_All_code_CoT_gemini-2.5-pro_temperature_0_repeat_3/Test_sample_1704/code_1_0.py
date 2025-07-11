import math
from fractions import Fraction

# A list to store all found solutions (sets of integers)
solutions = []

def find_solutions(n, target, start_x, path):
    """
    Recursively finds sets of n distinct integers whose reciprocals sum to target.
    
    Args:
        n (int): The number of integers remaining to find.
        target (Fraction): The target sum of reciprocals.
        start_x (int): The minimum value for the next integer to ensure distinctness and order.
        path (list): The list of integers found so far in the current path.
    """
    # Base case: we only need to find one more number.
    if n == 1:
        # The last number x must satisfy 1/x = target.
        # This means x = 1/target.
        # For x to be an integer, target's numerator must be 1.
        if target.numerator == 1:
            x = target.denominator
            # The number must also be larger than the previous one in the path.
            if x >= start_x:
                solution = path + [x]
                solutions.append(solution)
        return

    # Recursive step: find the next number in the sequence.
    # We establish a search range for the current integer 'x'.
    # Lower bound is start_x.
    # Upper bound: since all n remaining numbers are >= x, n/x >= target, so x <= n/target.
    upper_bound = math.floor(n / target)
    
    for x in range(start_x, upper_bound + 1):
        new_target = target - Fraction(1, x)
        
        # Pruning: if the new target is non-positive, any larger x will also result
        # in a non-positive target, so we can stop searching in this loop.
        if new_target <= 0:
            break
            
        # Recurse to find the remaining n-1 numbers.
        find_solutions(n - 1, new_target, x + 1, path + [x])

def solve_and_print():
    """
    Main function to solve the problem for n=4 and print the results.
    """
    N = 4
    # Start the search for 4 numbers whose reciprocals sum to 1.
    # The first number must be at least 2.
    find_solutions(N, Fraction(1, 1), 2, [])
    
    if not solutions:
        print("No solutions found for n=4.")
        return

    print(f"Found {len(solutions)} sets of {N} distinct positive integers whose reciprocals sum to 1.")
    print("The sums of these sets form the elements of S(4):")
    
    # Calculate the sum for each solution set.
    s4_elements = []
    # Sort solutions by their sum for a clear presentation.
    for s in sorted(solutions, key=sum):
        current_sum = sum(s)
        s4_elements.append(current_sum)
        print(f"  Set {s} -> Sum = {current_sum}")
    
    # T(4) is the sum of the elements in S(4).
    # Since a set contains unique elements, we take the sum of the unique sums.
    unique_s4_elements = sorted(list(set(s4_elements)))
    t4 = sum(unique_s4_elements)
    
    print("\nThe set S(4) is the set of these unique sums:")
    print(f"S(4) = {set(unique_s4_elements)}")
    
    print("\nT(4) is the sum of all elements in S(4).")
    # Format the final equation as requested.
    final_equation = " + ".join(map(str, unique_s4_elements))
    print(f"T(4) = {final_equation} = {t4}")
    
    # The final answer in the required format
    print(f"\n<<<{t4}>>>")

# Run the solver
solve_and_print()