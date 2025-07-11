import sys

def solve():
    """
    Finds the integer x that minimizes the total length and prints the result.
    The input is read from stdin, with each line containing one integer.
    """
    try:
        a_list = [int(line) for line in sys.stdin]
        if not a_list:
            print("No input provided.")
            return
    except (ValueError, IndexError):
        # Example data if stdin is empty or invalid.
        # This part can be adjusted to handle different input methods.
        print("Reading from stdin failed or empty. Using example data: a_list = [10, 20, 31]")
        a_list = [10, 20, 31]

    n = len(a_list)
    sum_a = sum(a_list)
    max_a = 0
    for val in a_list:
        if val > max_a:
            max_a = val

    if max_a == 0:
        print("Optimal x can be any positive integer, minimal length is 0.")
        return

    # O(n + A) to create counts
    counts = [0] * (max_a + 1)
    for val in a_list:
        counts[val] += 1

    # O(A) to create num_ge (number of elements >= k)
    num_ge = [0] * (max_a + 2)
    for i in range(max_a, -1, -1):
        num_ge[i] = num_ge[i + 1] + counts[i]

    # O(A log A) to calculate all S(v)
    # S[v] = sum(floor(a_i / v) for all i)
    S = [0] * (max_a + 2)
    for v in range(1, max_a + 1):
        # S(v) = sum_{k=1}^{floor(A/v)} num_ge[k*v]
        for j in range(v, max_a + 1, v):
            S[v] += num_ge[j]

    min_len = sum_a # This is the length for x=1 or x > max_a
    best_x = 1

    # O(A) to find the best x
    for x in range(2, max_a + 2):
        # length = sum_a + (1-x) * S(x)
        current_len = sum_a + (1 - x) * S[x]
        if current_len < min_len:
            min_len = current_len
            best_x = x

    # The equation for the minimal total length L is sum_{i=1 to n} (floor(a_i/x) + a_i % x)
    # Printing the two key numbers of the final result: the optimal x and the minimal length.
    print(f"Optimal x: {best_x}")
    print(f"Minimal total length: {min_len}")

solve()