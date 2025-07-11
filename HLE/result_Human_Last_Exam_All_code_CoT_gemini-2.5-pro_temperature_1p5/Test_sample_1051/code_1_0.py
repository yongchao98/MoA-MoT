import sys

def solve():
    """
    This function finds the positive integer x that minimizes the total length
    for a given list of integers.
    """
    lines = sys.stdin.readlines()
    if not lines:
        # Handle case with no input
        return

    try:
        a = [int(line.strip()) for line in lines]
    except ValueError:
        print("Error: All inputs must be integers.")
        return

    if not a:
        return
        
    n = len(a)
    max_a = 0
    sum_a = 0
    for val in a:
        if val > max_a:
            max_a = val
        sum_a += val

    # If all a_i are 0, any x>=1 gives length 0. x=1 is a valid answer.
    if max_a == 0:
        print(1)
        return

    # counts[k] = number of times k appears in the input 'a'
    counts = [0] * (max_a + 1)
    for val in a:
        counts[val] += 1

    # num_ge[k] = number of a_i's such that a_i >= k
    num_ge = [0] * (max_a + 2)
    for i in range(max_a, 0, -1):
        num_ge[i] = num_ge[i + 1] + counts[i]

    # S[x] = sum(floor(a_i / x)) for all i
    # We can compute this for all x from 1 to max_a efficiently.
    # S[x] = sum over j=1... of num_ge[j*x]
    s_values = [0] * (max_a + 1)
    for x in range(1, max_a + 1):
        total_quotient = 0
        for y in range(x, max_a + 1, x):
            total_quotient += num_ge[y]
        s_values[x] = total_quotient
        
    min_len = float('inf')
    best_x = 1

    # Iterate through all possible x and find the minimum total length
    for x in range(1, max_a + 1):
        # Total_Length(x) = sum(a_i) + (1-x) * S(x)
        current_len = sum_a + (1 - x) * s_values[x]
        
        if current_len < min_len:
            min_len = current_len
            best_x = x

    # According to the problem "Remember in the final code you still need to
    # output each number in the final equation!". It is ambiguous, but the most
    # likely interpretation is to output the value of x that was found.
    print(best_x)

solve()
