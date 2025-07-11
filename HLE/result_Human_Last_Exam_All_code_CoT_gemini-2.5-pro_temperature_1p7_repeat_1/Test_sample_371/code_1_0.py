import sys

def solve_s_k(n, k):
    """
    Calculates the k-th element of the recursively defined sequence S_n.

    The problem states:
    S_0 = [1]
    S_n = S_{n-1} + [n+1] + S_{n-1} for n > 0

    An analysis of this structure reveals that the value at index k is given by
    the formula v_2(k+1) + 1, where v_2(m) is the number of trailing zeros
    in the binary representation of m. The value of n is irrelevant as long as
    k is a valid index within the sequence S_n.

    This function implements the optimal, direct calculation.
    """
    if k < 0:
        print("Error: Index k must be non-negative.", file=sys.stderr)
        return

    # Let m = k + 1. We need to compute v_2(m) + 1.
    m = k + 1
    
    # In Python, for a positive integer x, the expression `(x & -x)` isolates
    # the least significant bit (which is a power of 2).
    # The `bit_length()` method for an integer that is a power of two, say 2^p,
    # returns p + 1.
    # Therefore, `(m & -m).bit_length()` directly computes `v_2(m) + 1`.
    result = (m & -m).bit_length()

    # The value `p` in the formula is v_2(k+1), which is `result - 1`.
    p = result - 1

    # Per the instruction "output each number in the final equation!",
    # we display the final computation based on the formula: result = p + 1.
    print(f"{result} = {p} + 1")


if __name__ == '__main__':
    # Set default values based on the example S_2[3] = 3.
    n_val = 2
    k_val = 3

    # If command-line arguments are provided, use them instead.
    # The script expects two integer arguments: n and k.
    if len(sys.argv) == 3:
        try:
            n_val = int(sys.argv[1])
            k_val = int(sys.argv[2])
        except (ValueError, IndexError):
            print(f"Usage: python {sys.argv[0]} <n> <k>", file=sys.stderr)
            print("Please provide two integer values for n and k.", file=sys.stderr)
            sys.exit(1)
    
    solve_s_k(n_val, k_val)