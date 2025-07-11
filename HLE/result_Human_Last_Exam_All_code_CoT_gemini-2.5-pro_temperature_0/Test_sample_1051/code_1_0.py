import sys
import math

def solve():
    """
    Reads integers from standard input, finds an integer x that minimizes the
    total length, and prints the final equation.
    """
    try:
        lines = sys.stdin.readlines()
        if not lines:
            return
        a = [int(line.strip()) for line in lines]
        n = len(a)
        if n == 0:
            return
    except (IOError, ValueError) as e:
        # In case of input errors, you might want to handle them gracefully.
        # For this problem, we assume valid integer inputs.
        return

    A = 0
    for val in a:
        if val > A:
            A = val

    if A == 0:
        # If all a_i are 0, any x > 0 gives length 0. The smallest positive x is 1.
        final_parts = [f"floor(0 / 1) + (0 % 1)" for _ in range(n)]
        equation = " + ".join(final_parts)
        print(f"{equation} = 0")
        return

    sum_a = sum(a)
    
    # diff_S[x] will store S(x) - S(x-1), where S(x) = sum(floor(a_i / x))
    diff_S = [0] * (A + 2)

    for val in a:
        if val == 0:
            continue
        
        # The contribution of `val` to S(x) is floor(val / x).
        # We add its difference function, d(x) = floor(val/x) - floor(val/(x-1)),
        # to the total diff_S array.
        
        # At x=1, the value is `val`. The change from S(0) (assumed 0) is `val`.
        diff_S[1] += val
        
        # The value floor(val/x) only changes at specific breakpoints.
        # We find these breakpoints and update the diff_S array accordingly.
        limit = int(math.sqrt(val))
        breakpoints = set()
        for i in range(1, limit + 1):
            # Breakpoints occur at x and floor(val/i) + 1
            breakpoints.add(i + 1)
            breakpoints.add(val // i + 1)
        
        for x in sorted(list(breakpoints)):
            if x > val + 1:
                continue
            # Calculate the change in floor(val/x) at this breakpoint
            change = (val // x) - (val // (x - 1))
            if change != 0:
                diff_S[x] += change

    min_len = float('inf')
    best_x = -1
    
    current_S = 0
    # Iterate from x=1 to A+1 to find the optimal x
    for x in range(1, A + 2):
        current_S += diff_S[x]
        
        # Total length L(x) = sum(a_i) - (x-1) * S(x)
        current_len = sum_a - (x - 1) * current_S
        
        if current_len < min_len:
            min_len = current_len
            best_x = x

    # Format and print the final equation as requested
    final_parts = []
    for val in a:
        final_parts.append(f"floor({val} / {best_x}) + ({val} % {best_x})")
    
    equation = " + ".join(final_parts)
    
    # Recalculate total_len as an integer for the final output
    total_len = 0
    for val in a:
        total_len += (val // best_x) + (val % best_x)
        
    print(f"{equation} = {total_len}")

solve()