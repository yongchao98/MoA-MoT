import sys

def solve():
    """
    Finds the positive integer x that minimizes the total length for a list of integers a_i.
    The length for each a_i is defined as floor(a_i/x) + (a_i mod x).
    """
    # Step 1: Read input integers and find the maximum value A.
    try:
        lines = sys.stdin.readlines()
        if not lines:
            print("No input provided. Please provide integers, one per line.")
            return
        a = [int(line.strip()) for line in lines]
    except (ValueError, IndexError):
        print("Invalid input. Please provide one integer per line.")
        return

    if not a:
        print("Input list is empty.")
        return

    A = 0
    for val in a:
        if val > A:
            A = val

    # If all a_i are 0, any x > 0 gives length 0. x=1 is a valid choice.
    if A == 0:
        print("The optimal value for x is 1")
        print("The minimized total length is calculated as: 0 = 0")
        return

    # Step 3: Create a frequency array `counts`.
    counts = [0] * (A + 1)
    for val in a:
        counts[val] += 1

    # Step 4: Create a suffix sum array `num_ge`.
    num_ge = [0] * (A + 2)
    for v in range(A, -1, -1):
        num_ge[v] = num_ge[v+1] + counts[v]

    # Step 5-9: Find the best x that maximizes g(x) = (x-1) * S(x).
    # g(1) = 0.
    max_g = 0
    best_x = 1

    # Iterate x from 2 to A. For x > A, g(x) = 0.
    for x in range(2, A + 1):
        # Calculate S(x) = sum(floor(a_i / x)) using the pre-computed array.
        s_x = 0
        for kx in range(x, A + 1, x):
            s_x += num_ge[kx]
        
        # Calculate g(x) = (x-1) * S(x)
        g_x = (x - 1) * s_x
        
        if g_x > max_g:
            max_g = g_x
            best_x = x

    # Step 11: Calculate the minimal total length with the best_x found and show the equation.
    print(f"The optimal value for x is {best_x}")
    
    min_total_length = 0
    length_parts = []
    for val in a:
        q = val // best_x
        r = val % best_x
        length_i = q + r
        min_total_length += length_i
        length_parts.append(f"({q} + {r})")

    equation_str = " + ".join(length_parts)
    print(f"The minimized total length is calculated as: {equation_str} = {min_total_length}")

if __name__ == "__main__":
    solve()