import math

def get_kth_element(n: int, k: int) -> None:
    """
    Calculates the k-th element of a recursively defined sequence S_n.

    The sequence is defined as:
    S_0 = [1]
    S_n = S_{n-1} + [n+1] + S_{n-1} for n > 0

    The optimal solution is derived from the observation that the k-th element
    of this sequence (for a sufficiently large n) is determined by the
    number of trailing zeros in the binary representation of k+1.
    This value is also known as the 2-adic valuation, v_2(k+1).
    The formula is: S_n[k] = v_2(k+1) + 1.

    This function implements this direct formula, which is highly efficient
    and handles the large integer constraints specified in the problem.

    Args:
        n: The sequence number. Assumed to be large.
        k: The 0-based index of the element to find. Assumed to be large.
    """
    if k < 0:
        print("Error: Index k must be non-negative.")
        return

    # Let x = k+1. We need to calculate v_2(x), the number of trailing zeros in x's binary form.
    x = k + 1

    # The bitwise hack `x & -x` isolates the least significant 1-bit of x.
    # The result is a power of 2, say 2^p, where p = v_2(x).
    lsb_power_of_2 = x & -x

    # For a number y that is a power of two (y = 2^p), its bit length is p+1.
    # So, we can find p by computing `y.bit_length() - 1`.
    p = lsb_power_of_2.bit_length() - 1

    # The value of the element is v_2(k+1) + 1, which is p + 1.
    result = p + 1

    # The problem asks to output the numbers in the final equation.
    # The equation is S_n[k] = v_2(k+1) + 1.
    print(f"S_{n}[{k}] = v_2({k+1}) + 1 = {p} + 1 = {result}")

# Example from the problem description:
# get_kth_element(n=2, k=3)

# Example with large numbers
# This demonstrates that the function works with numbers larger than a 64-bit integer.
# try:
#     get_kth_element(n=10**20, k=(10**20) - 8)
# except TypeError:
#     # In some environments, int may not be usable for such large powers directly
#     # but the function itself supports large integer inputs.
#     pass