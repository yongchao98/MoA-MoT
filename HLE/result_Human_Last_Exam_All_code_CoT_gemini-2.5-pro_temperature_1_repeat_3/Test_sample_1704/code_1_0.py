import math
from fractions import Fraction

def find_t_n(n):
    """
    Finds T(n), the sum of all elements in S(n).
    S(n) is the set of sums of n distinct positive integers whose reciprocals sum to 1.
    """
    solutions = []

    def find_reciprocal_sums(target, k, min_val, path):
        """
        Recursively finds sets of k distinct integers > min_val whose reciprocals sum to target.
        - target: The Fraction object representing the sum to be achieved.
        - k: The number of integers left to find.
        - min_val: The minimum value for the next integer (must be > the previous one).
        - path: The list of integers found so far.
        """
        # Base case: we need to find the last term.
        if k == 1:
            # The reciprocal of the target must be an integer.
            if target.numerator == 1 and target.denominator > min_val:
                last_term = target.denominator
                solutions.append(path + [last_term])
            return

        # Determine the search range for the next term, i.
        # Lower bound for i:
        # i must be greater than min_val.
        # Also, target > 1/i => i > 1/target.
        low = max(min_val + 1, math.floor(1 / target) + 1)

        # Upper bound for i:
        # target < k/i => i < k/target.
        high = math.floor(k / target)

        for i in range(low, high + 1):
            new_target = target - Fraction(1, i)
            # If new_target is non-positive, further search is futile.
            if new_target > 0:
                find_reciprocal_sums(new_target, k - 1, i, path + [i])

    # Initial call for n=4, target=1. The first number must be at least 2.
    # We set min_val=1, so the first number starts at 1+1=2.
    find_reciprocal_sums(Fraction(1, 1), n, 1, [])
    
    # S(n) contains the sum of each solution set. Using a set to handle duplicates if any.
    s_n = {sum(sol) for sol in solutions}
    
    # T(n) is the sum of all elements in S(n).
    t_n = sum(s_n)
    
    return solutions, sorted(list(s_n)), t_n

def main():
    """
    Main function to solve the problem for n=4 and print the results.
    """
    n = 4
    solutions, sorted_sums, t_4 = find_t_n(n)
    
    # Sort solutions based on their sum for a cleaner presentation.
    sorted_solutions = sorted(solutions, key=sum)

    print(f"To find T(4), we first identify all sets of {n} distinct positive integers whose reciprocals sum to 1.")
    print("The discovered sets {x1, x2, x3, x4} are:")
    for sol in sorted_solutions:
        equation = " + ".join([f"1/{x}" for x in sol])
        print(f"- {{{', '.join(map(str, sol))}}}: Sum = {sum(sol):<2}  (since {equation} = 1)")

    print(f"\nThe set of sums, S(4), is {{{', '.join(map(str, sorted_sums))}}}.")
    
    calculation_str = " + ".join(map(str, sorted_sums))
    print(f"\nT(4) is the sum of the elements in S(4).")
    print(f"T(4) = {calculation_str} = {t_4}")

if __name__ == "__main__":
    main()