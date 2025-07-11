import sys

def solve():
    """
    Finds the optimal positive integer x that minimizes the total length.
    The total length is the sum over all input integers a_i of:
    floor(a_i / x) + (a_i mod x)
    """
    try:
        # For local testing, you can provide a list directly.
        # For example: a = [6, 10, 20]
        # In a typical competitive programming setup, we read from stdin.
        lines = sys.stdin.readlines()
        if not lines:
            # Handle empty input
            a = []
        else:
            # Assumes the first line is n and the rest are a_i, or all lines are a_i.
            # This parsing handles inputs where numbers are on separate lines.
            a = [int(line.strip()) for line in lines if line.strip().isdigit()]

    except (IOError, ValueError):
        # Fallback for environments where stdin is not available or if parsing fails.
        # You can hardcode an example here.
        a = [6, 10, 20]

    if not a:
        print("Input list is empty. Optimal x is 1, total length is 0.")
        return
        
    n = len(a)
    A = 0
    s_a = 0
    for val in a:
        if val > A:
            A = val
        s_a += val

    # O(n + A)
    counts = [0] * (A + 1)
    for val in a:
        counts[val] += 1

    # O(A)
    num_ge = [0] * (A + 2)
    for i in range(A, 0, -1):
        num_ge[i] = num_ge[i+1] + counts[i]

    # O(A log A)
    s_q = [0] * (A + 2)  # s_q[x] will store sum(floor(a_i/x))
    for x in range(1, A + 1):
        val = 0
        # This computes s_q(x) in O(A/x) time
        for k in range(1, A // x + 1):
            val += num_ge[k * x]
        s_q[x] = val

    best_x = 1
    # For x=1, length is s_a + (1-1)*s_q[1] = s_a
    min_len = float(s_a)

    # O(A)
    for x in range(2, A + 1):
        current_len = s_a + (1 - x) * s_q[x]
        if current_len < min_len:
            min_len = current_len
            best_x = x

    # Final Output
    x = best_x
    print(f"The optimal value for x is {x}.")
    print("The minimized total length is calculated as follows:")
    
    equation_parts = []
    final_total_len = 0
    for val in a:
        q = val // x
        r = val % x
        length_i = q + r
        final_total_len += length_i
        equation_parts.append(f"({q} + {r})")
    
    print(" + ".join(equation_parts) + f" = {int(final_total_len)}")

solve()
