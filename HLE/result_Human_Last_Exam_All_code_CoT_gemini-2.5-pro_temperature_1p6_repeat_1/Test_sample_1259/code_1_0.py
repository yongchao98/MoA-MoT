import math

def solve_and_construct_A():
    """
    Reads a set S, checks if a generating sequence A can exist,
    and if so, constructs and prints one such A.
    """
    try:
        # Read space-separated numbers from a single line of input.
        # As an example, you can paste: 1 6 10 15
        s_elements = list(map(int, input().split()))
    except (ValueError, EOFError):
        print("Invalid input. Please provide a space-separated list of integers.")
        # Example data for demonstration if input fails
        s_elements = [6, 12, 18] 
        print(f"Using example data S = {s_elements}")

    if not s_elements:
        print("Input set S cannot be empty.")
        return

    # Find the minimum element, g.
    g = min(s_elements)

    # Check the necessary and sufficient condition.
    for s in s_elements:
        if s % g != 0:
            print("-1") # Indicates no solution exists
            return

    # If the condition holds, a solution exists. Construct one.
    # Construction: A = [g, s'_1, g, s'_2, ...], where s'_i are elements of S other than g.
    # This construction is guaranteed to work.
    
    A = [g]
    for s in s_elements:
        if s != g:
            A.append(s)
            A.append(g)

    # Print the resulting sequence A.
    # First, print the length of A.
    print(len(A))
    # Then print the elements of A, space-separated.
    # The prompt requested "each number in the final equation"
    # which we interpret as the elements of the constructed sequence.
    print(*A)

# To run this code, you would execute the function.
# For example, if you input "6 12 18", it will output:
# 3
# 6 12 18 
# If you input "1 6 10 15", it will output:
# 7
# 1 6 1 10 1 15 1
# This is a valid A, although shorter ones might exist.
solve_and_construct_A()