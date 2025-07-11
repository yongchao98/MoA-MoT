import sys

def solve():
    """
    Finds the positive integer x that minimizes the total length and prints the result.
    The total length for a list of integers a_i is Sum(floor(a_i/x) + (a_i mod x)).
    """
    # Step 1: Read input from stdin
    try:
        lines = sys.stdin.readlines()
        if not lines:
            print("No input provided.")
            return
        a = [int(line.strip()) for line in lines]
        n = len(a)
    except (ValueError, IndexError):
        print("Invalid input. Please provide one integer per line.")
        return

    if n == 0:
        print("Input is empty.")
        return

    # Step 2: Find the maximum value A in the input list a
    A = 0
    for val in a:
        if val > A:
            A = val

    if A == 0:
        # If all a_i are 0, any x > 0 gives a total length of 0. We can pick x=1.
        final_lengths = [0] * n
        total_length = 0
        print(" + ".join(map(str, final_lengths)) + " = " + str(total_length))
        return

    # Step 3: Count occurrences of each number
    counts = [0] * (A + 1)
    for val in a:
        counts[val] += 1

    # Step 4: Compute suffix sums S, where S[v] = count of numbers >= v
    S = [0] * (A + 2)
    for v in range(A, -1, -1):
        S[v] = S[v + 1] + counts[v]

    # Step 5: Iterate through possible x values to find the one that maximizes F(x)
    best_x = 1
    max_F = 0  # F(1) is 0

    for x in range(2, A + 1):
        # Calculate Q(x) = Sum(floor(a_i / x)) efficiently using the S array
        Q_x = 0
        for k in range(1, A // x + 1):
            Q_x += S[k * x]
        
        # Calculate F(x) = (x-1) * Q(x)
        F_x = (x - 1) * Q_x
        
        if F_x > max_F:
            max_F = F_x
            best_x = x

    # Step 6: Calculate the final lengths and total length with the best x
    final_lengths = []
    total_length = 0
    for val in a:
        length_i = (val // best_x) + (val % best_x)
        final_lengths.append(length_i)
        total_length += length_i
    
    # Print the result in the format "l_1 + l_2 + ... = total"
    print(" + ".join(map(str, final_lengths)) + " = " + str(total_length))

if __name__ == '__main__':
    solve()