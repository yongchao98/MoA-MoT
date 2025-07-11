import sys

def solve():
    """
    Finds the positive integer x that minimizes the total length for a given list of integers.
    """
    # Step 1: Read all integers a_i from standard input
    a_values = []
    try:
        for line in sys.stdin:
            line = line.strip()
            if line:
                val = int(line)
                if val > 0:
                    a_values.append(val)
    except ValueError:
        print("Invalid input. Please provide only integers.", file=sys.stderr)
        return

    if not a_values:
        print("No positive integers provided.")
        return

    # Step 2: Determine maximum value A and sum
    A = max(a_values)
    sum_a = sum(a_values)

    # Step 3: Efficiently compute Q(x) for x in [1, A]
    # Create counts array
    counts = [0] * (A + 1)
    for val in a_values:
        counts[val] += 1
    
    # Create num_ge array (number of a_i >= v) using a suffix sum on counts
    num_ge = [0] * (A + 2)
    for v in range(A, 0, -1):
        num_ge[v] = num_ge[v + 1] + counts[v]
        
    # Compute Q[x] for all x from 1 to A
    Q = [0] * (A + 1)
    for x in range(1, A + 1):
        q_val = 0
        for multiple_x in range(x, A + 1, x):
            q_val += num_ge[multiple_x]
        Q[x] = q_val
    
    # Step 4: Find the best x that minimizes the total length
    min_len = sum_a  # This corresponds to the length for x=1
    best_x = 1

    for x in range(2, A + 1):
        current_len = sum_a + (1 - x) * Q[x]
        if current_len < min_len:
            min_len = current_len
            best_x = x

    # Step 5: Print the results with the optimal x
    print(f"The optimal value for x is {best_x}.")
    print("This minimizes the total length, with the following breakdown:")
    
    total_length = 0
    for i, a in enumerate(a_values, 1):
        quotient = a // best_x
        remainder = a % best_x
        length_i = quotient + remainder
        total_length += length_i
        print(f"length_{i} = floor({a}/{best_x}) + ({a} % {best_x}) = {quotient} + {remainder} = {length_i}")
    
    print(f"\nMinimal total length = {int(total_length)}")

if __name__ == '__main__':
    solve()