import math
import sys

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of integers."""
    if not numbers:
        return 0
    # All numbers are positive, so result will be positive.
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = math.gcd(result, numbers[i])
    return result

def solve():
    """
    Reads a set S from standard input.
    If a generating sequence A exists, it prints one such valid sequence.
    Otherwise, it prints -1.
    """
    try:
        # Read the set S from a single line of space-separated integers.
        s_input = sys.stdin.readline()
        if not s_input.strip():
            print("-1")
            return
        
        # Create a set to handle duplicates, then convert to a list.
        S_set = set(map(int, s_input.split()))
        S_list = list(S_set)
        
        if not S_list or any(x <= 0 for x in S_list):
            print("-1")
            return

    except (ValueError, IndexError):
        # Handle cases where input is not a list of integers.
        print("-1")
        return

    # A solution exists if and only if min(S) == gcd(S).
    g = gcd_list(S_list)
    min_s = min(S_list)

    if g != min_s:
        # If the condition is not met, no solution exists.
        print("-1")
        return

    # If the condition is met, a solution is guaranteed. We construct one.
    # S_prime contains all elements of S except for g, sorted.
    S_prime = sorted([s for s in S_set if s != g])

    # The constructed sequence A.
    A = []
    if not S_prime:
        # This case occurs if S only contains one element, g.
        A = [g]
    else:
        # A robust construction is A = [s_1, g, s_2, g, ..., g, s_k].
        # This prevents generating unwanted GCDs.
        for i, s_val in enumerate(S_prime):
            A.append(s_val)
            if i < len(S_prime) - 1:
                A.append(g)

    # Print the elements of the constructed sequence A, separated by spaces.
    print(*A)

# To run this code, you would typically call the solve() function.
# For example:
# if __name__ == "__main__":
#     solve()
# Input (e.g., 1 6 10 15) would be read from stdin.
# Output for that input would be: 6 1 10 1 15