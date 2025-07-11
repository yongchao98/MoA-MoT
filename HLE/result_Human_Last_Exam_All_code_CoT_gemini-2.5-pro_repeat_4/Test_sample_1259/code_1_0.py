import math
import sys

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of numbers."""
    if not numbers:
        return 0
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = math.gcd(result, numbers[i])
    return result

def solve():
    """
    Reads a set S, checks if a solution A exists, and if so,
    constructs and prints one.
    """
    try:
        # Read the input set S from a single line of space-separated integers
        s_str = input().strip().split()
        if not s_str:
            print(-1)
            return
        S = sorted([int(x) for x in s_str])
    except (IOError, ValueError):
        # Handle potential empty input or non-integer values
        print("Invalid input. Please provide a list of positive integers.")
        # Example data for demonstration if input fails
        # S = [6, 12, 30] 
        # S = {1, 6, 10, 15}
        return

    # The minimum element of S
    m = S[0]
    
    # The gcd of all elements in S
    g = gcd_list(S)

    # Check the necessary and sufficient condition
    if m != g:
        print(-1)
        return

    # If the condition holds, a solution is guaranteed to exist.
    # We construct one by interleaving the minimum element 'm'
    # with the other elements of S.
    A = []
    if len(S) == 1:
        A = S
    else:
        # Construction: [s'_1, m, s'_2, m, ..., s'_k]
        # This has length 2 * |S'| - 1 = 2 * (|S|-1) - 1 = 2|S| - 3
        s_prime = S[1:]
        for i, val in enumerate(s_prime):
            A.append(val)
            # Don't add 'm' after the last element
            if i < len(s_prime) - 1:
                A.append(m)
    
    # Another simple construction is to interleave all elements with m
    # which gives a slightly longer but still valid sequence.
    # Let's use that for simplicity of implementation.
    # A = [s'_1, m, s'_2, m, ... s'_k, m]
    A_interleaved = []
    if len(S) > 1:
        for val in S[1:]:
            A_interleaved.append(val)
            A_interleaved.append(S[0])
    else:
        A_interleaved = S
        
    # The problem asks to output the result, not ask the user to copy.
    # The prompt also mentions "final equation", which is ambiguous.
    # A standard output format for sequences is length followed by elements.
    print(len(A_interleaved))
    # The * operator unpacks the list elements for the print function
    print(*A_interleaved)

# To run the solution, uncomment the line below or run this script in a
# terminal and provide input.
# solve()

# Example usage for demonstration:
# For S = {2, 6, 10}, a valid A is [6, 2, 10]
# For S = {1, 6, 10, 15}, a valid A is [6, 1, 10, 1, 15]
# Let's use the code's construction for S = {6, 12, 18}
# S (sorted) = [6, 12, 18]
# m = 6
# S' = [12, 18]
# A = [12, 6, 18]
# The code outputs:
# 4
# 12 6 18 6
# which is also a valid construction.
# Let's run with the code that generates A=[s'1,m,s'2,m...s'k,m]
# input: 6 12 18
# S = [6, 12, 18]. m = 6. S' = [12, 18]
# A = [] -> A.append(12) -> A.append(6) -> A.append(18) -> A.append(6)
# A = [12, 6, 18, 6].
# len(A) = 4.
# Output:
# 4
# 12 6 18 6
# This is a valid sequence.
# gcd(12)=12, gcd(6)=6, gcd(18)=18, gcd(12,6)=6, gcd(6,18)=6, gcd(12,6,18)=6
# S_A = {6, 12, 18}. It works.
