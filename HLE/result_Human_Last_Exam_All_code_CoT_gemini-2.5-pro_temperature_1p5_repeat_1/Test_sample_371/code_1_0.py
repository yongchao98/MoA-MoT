def solve_kth_element():
    """
    Calculates the k-th element of S_n using the optimal method.

    The problem specifies that n and k are very large (`> INT_MAX / 2`), so we
    use numbers that require arbitrary-precision arithmetic, which Python handles
    natively. The chosen formula S_n[k] = ctz(k+1) + 1 is independent of n.
    """
    # Example large values for n and k that satisfy the problem's constraints.
    # We choose k = 10**100 - 1 to get a non-zero result for ctz.
    # 10**100 = (2*5)**100 = 2**100 * 5**100, so ctz(10**100) is 100.
    n = 10**100
    k = 10**100 - 1

    # The value to analyze is k + 1
    val = k + 1

    # The count of trailing zeros (ctz) is found by isolating the
    # least significant bit (LSB) and finding its position.
    # (val & -val) isolates the LSB.
    # .bit_length() - 1 gives the log base 2 of the LSB, which is the ctz.
    lsb = val & -val
    num_trailing_zeros = lsb.bit_length() - 1
    
    # The final result is ctz + 1
    result = num_trailing_zeros + 1
    
    # The prompt requires outputting each number in the final equation.
    # The final calculation is `num_trailing_zeros + 1 = result`.
    print(f"{num_trailing_zeros} + 1 = {result}")

# Execute the function to solve the task.
solve_kth_element()