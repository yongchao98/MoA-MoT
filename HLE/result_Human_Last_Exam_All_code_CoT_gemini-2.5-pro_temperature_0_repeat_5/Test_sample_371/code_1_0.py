def solve_kth_element():
    """
    Calculates the k-th element of the recursively defined sequence S_n.

    The function uses a direct mathematical formula derived from the sequence's properties,
    making it highly efficient and independent of the large value of n.
    """
    # Example values as given in the prompt.
    # The solution is independent of n, as long as k is a valid index.
    # The large value of n given in the problem description serves to invalidate
    # naive solutions that would iterate or recurse up to n.
    n = 2
    k = 3

    # The value of S_n[k] is equal to the number of trailing zeros
    # in the binary representation of (k+1), plus one.
    # Let x = k + 1.
    # The number of trailing zeros can be found using a bitwise trick.
    # `x & -x` isolates the lowest set bit of x. Let this be `lsb`.
    # For example, if x=12 (1100), lsb=4 (0100).
    # `lsb` will always be a power of 2. If lsb = 2^p, then p is the
    # number of trailing zeros.
    # In Python, for a power of two `p`, `p.bit_length()` returns log2(p) + 1.
    # So, the number of trailing zeros is `lsb.bit_length() - 1`.
    # The value of the element is ctz + 1, which simplifies to `lsb.bit_length()`.

    # We handle the case k=-1 (x=0) for robustness, though problem constraints
    # ensure k is large and positive.
    x = k + 1
    if x <= 0:
        # This case is not expected based on problem description (k > 0)
        result = "Invalid k"
    else:
        lowest_set_bit = x & -x
        result = lowest_set_bit.bit_length()

    # The prompt requires printing the numbers in the final equation S_n[k] = result
    print(f"S_{n}[{k}] = {result}")


if __name__ == '__main__':
    solve_kth_element()