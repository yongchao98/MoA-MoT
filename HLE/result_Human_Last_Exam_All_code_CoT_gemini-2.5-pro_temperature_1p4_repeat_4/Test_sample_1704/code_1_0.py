import math
from fractions import Fraction

def solve_T4():
    """
    Finds T(4), the sum of all elements in S(4).
    S(4) is the set of sums of 4 distinct positive integers whose reciprocals sum to 1.
    """
    
    # This list will store all found sets of integers (solutions).
    solutions = []

    def find_partitions(target_sum, k, min_denom, path):
        """
        Recursively finds sets of k distinct integers whose reciprocals sum to target_sum.
        
        Args:
            target_sum: The remaining sum to be achieved (as a Fraction).
            k: The number of integers left to find.
            min_denom: The minimum value for the next integer (to ensure distinctness).
            path: The list of integers found so far.
        """
        # Base case: We are looking for the last integer.
        if k == 1:
            # The reciprocal of the last integer must be equal to the remaining target_sum.
            # So, the integer itself is 1 / target_sum.
            if target_sum.numerator == 1:
                last_denom = target_sum.denominator
                # It must be a distinct integer, larger than the previous one.
                if last_denom > min_denom:
                    solutions.append(path + [last_denom])
            return

        # Recursive step: Determine the search range for the next integer 'd'.
        # Lower bound for d: It must be greater than the previous one.
        # It must also be large enough such that 1/d < target_sum.
        # d > 1/target_sum
        lower_bound = min_denom + 1
        if target_sum > 0:
             # ceil(1/target_sum)
            lower_bound = max(lower_bound, math.ceil(1 / target_sum))

        # Upper bound for d: It must be small enough such that k * (1/d) > target_sum.
        # d < k/target_sum
        upper_bound = math.floor(k / target_sum)

        for d in range(lower_bound, upper_bound + 1):
            new_target_sum = target_sum - Fraction(1, d)
            find_partitions(new_target_sum, k - 1, d, path + [d])

    # Initial call to start the search for n=4.
    # We are looking for 4 terms that sum to 1. The smallest possible integer is 2,
    # since if x1=1, the sum is already 1, which is impossible with other positive integers.
    # So the first integer must be > 1. We set min_denom=1, so the first loop starts at 2.
    find_partitions(Fraction(1, 1), 4, 1, [])

    print("Found all sets of 4 distinct positive integers whose reciprocals sum to 1:")
    
    # S(4) is a set, so we use a set to store the unique sums.
    s4 = set()
    
    # Sort solutions for consistent output
    solutions.sort()

    for sol in solutions:
        current_sum = sum(sol)
        s4.add(current_sum)
        
        reciprocal_eq = " + ".join([f"1/{x}" for x in sol])
        sum_eq = " + ".join([str(x) for x in sol])
        
        print(f"{reciprocal_eq} = 1, sum = {sum_eq} = {current_sum}")
    
    # T(4) is the sum of all elements in S(4).
    t4 = sum(s4)
    
    print("\nThe set of sums, S(4), is:", sorted(list(s4)))
    print("The sum of these elements, T(4), is:", t4)
    return t4

if __name__ == '__main__':
    final_answer = solve_T4()
    # The final answer in the requested format
    print(f"\n<<<{final_answer}>>>")
