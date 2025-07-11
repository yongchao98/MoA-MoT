import math
from fractions import Fraction

def find_t_n(n_val):
    """
    Finds the value of T(n), the sum of all elements in S(n).
    S(n) is the set of sums of n distinct positive integers whose reciprocals sum to 1.
    """

    # This set will store the unique tuples of solutions
    solutions = set()

    def find_reciprocal_sums(n, target_sum, min_x, terms):
        """
        Recursively finds sets of n distinct integers whose reciprocals sum to a target.

        Args:
            n (int): The number of distinct integers still to find.
            target_sum (Fraction): The target sum for the reciprocals of the remaining integers.
            min_x (int): The minimum value for the next integer to ensure they are distinct
                         and in ascending order.
            terms (list): The list of integers found so far in the current path.
        """
        # Base case: we need to find the last integer
        if n == 1:
            # The remaining target sum must be a unit fraction 1/x for the solution to be valid.
            if target_sum.numerator == 1:
                x = target_sum.denominator
                # The new integer must be larger than the previous one.
                if x > min_x:
                    # Found a valid solution set. Add it to our collection.
                    solution_set = tuple(terms + [x])
                    solutions.add(solution_set)
            return

        # --- Recursive step: find bounds for the next integer 'x' and search ---
        
        # Lower bound for x:
        # 1. x must be > than the previously found integer (min_x) to be distinct.
        # 2. 1/x must be < target_sum, so x > 1/target_sum.
        lower_bound = max(min_x + 1, math.floor(1 / target_sum) + 1)
        
        # Upper bound for x:
        # Since x is the smallest of the remaining n integers, target_sum < n * (1/x).
        # This implies x < n / target_sum.
        # The range function's stop is exclusive, so we use math.ceil for the upper limit.
        upper_bound = math.ceil(n / target_sum)

        for x in range(lower_bound, upper_bound):
            new_target_sum = target_sum - Fraction(1, x)
            # Only proceed if the remaining sum is positive.
            if new_target_sum > 0:
                find_reciprocal_sums(n - 1, new_target_sum, x, terms + [x])

    # Start the search for solutions to sum(1/x_i) = 1 for n_val integers.
    # The smallest possible x_i is 2, since if x_1=1, the sum is already >= 1.
    # So we can set the initial min_x to 1.
    find_reciprocal_sums(n_val, Fraction(1, 1), 1, [])

    # S(n) is the set of sums of the integer sets found.
    # Using a set automatically handles any duplicate sums.
    s_n = {sum(s) for s in solutions}

    # T(n) is the sum of all elements in S(n).
    t_n = sum(s_n)
    
    # --- Output the results as requested ---

    # Sort the unique sums for clear presentation.
    sorted_sums = sorted(list(s_n))
    equation_str = " + ".join(map(str, sorted_sums))

    print(f"Found {len(solutions)} sets of {n_val} distinct positive integers whose reciprocals sum to 1.")
    print(f"The unique sums of these sets are: {', '.join(map(str, sorted_sums))}")
    print("\nT(4) is the sum of these unique values:")
    print(f"{equation_str} = {t_n}")

if __name__ == '__main__':
    # Find T(n) for n=4
    find_t_n(4)