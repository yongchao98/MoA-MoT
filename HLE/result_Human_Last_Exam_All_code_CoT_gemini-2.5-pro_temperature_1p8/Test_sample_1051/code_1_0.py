import sys

def solve():
    """
    Reads integers, finds the optimal x to minimize total length,
    and prints the detailed calculation.
    """
    # Assuming input integers are provided one per line.
    # For local testing, you can use a hardcoded list:
    # a = [6, 10, 20]
    
    # To read from standard input:
    try:
        lines = sys.stdin.readlines()
        if not lines:
            print("Input is empty.")
            return
        a = [int(line.strip()) for line in lines]
    except (ValueError, IndexError):
        print("Invalid input. Please provide a list of integers, one per line.")
        return

    n = len(a)
    if n == 0:
        print("Input list is empty.")
        return

    max_a = 0
    sum_a = 0
    for val in a:
        if val > max_a:
            max_a = val
        sum_a += val

    if max_a == 0:
        parts = ["(0 + 0)"] * n
        print(f"{' + '.join(parts)} = 0")
        return
        
    counts = [0] * (max_a + 1)
    for val in a:
        counts[val] += 1
    
    suffix_counts = [0] * (max_a + 2)
    for i in range(max_a, 0, -1):
        suffix_counts[i] = suffix_counts[i+1] + counts[i]

    min_len = sum_a  # Length for x=1
    best_x = 1

    # We start from x=2 because x=1 is the initial case
    for x in range(2, max_a + 1):
        s_x = 0
        for j in range(x, max_a + 1, x):
            s_x += suffix_counts[j]
        
        current_len = sum_a - (x - 1) * s_x
        
        if current_len < min_len:
            min_len = current_len
            best_x = x

    # Output the final equation with the best x
    parts = []
    final_len = 0
    for val in a:
        q = val // best_x
        r = val % best_x
        parts.append(f"({q} + {r})")
        final_len += q + r

    print(f"{' + '.join(parts)} = {final_len}")

solve()