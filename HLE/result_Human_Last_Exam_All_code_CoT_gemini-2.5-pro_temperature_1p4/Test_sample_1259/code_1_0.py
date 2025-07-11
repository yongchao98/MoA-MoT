import math
from functools import reduce

def restore_A_from_S(s_list):
    """
    Given a set S, restores a possible original sequence A.

    Args:
        s_list: A list of integers representing the set S.

    Returns:
        A list representing the sequence A, or None if no solution exists.
    """
    if not s_list:
        return None

    S = sorted(list(set(s_list)))

    def gcd_list(numbers):
        if len(numbers) == 0:
            return 0
        return reduce(math.gcd, numbers)

    # Check the necessary and sufficient condition: min(S) == gcd(S)
    g = S[0]
    if g != gcd_list(S):
        return None

    # If the condition holds, a solution is guaranteed to exist.
    # We construct one such solution. A simple one is to interleave
    # the smallest element 'g' with the other elements of S.
    # For S = {g, s2, s3, ..., sk}, one valid A is [g, s2, g, s3, ..., g, sk].
    
    if len(S) == 1:
        return S

    A = []
    rest_of_S = S[1:]
    for val in rest_of_S:
        A.append(g)
        A.append(val)
        
    return A

if __name__ == '__main__':
    # Example usage:
    # Read input for S. Assuming space-separated integers on one line.
    try:
        s_input_str = input("Enter the set S as space-separated integers: ")
        s_input_list = list(map(int, s_input_str.split()))
    except (ValueError, IndexError):
        print("Invalid input. Please provide space-separated integers.")
        s_input_list = []

    if s_input_list:
        result_A = restore_A_from_S(s_input_list)
    
        if result_A is None:
            print(-1)
        else:
            # The problem asks to "output each number in the final equation!".
            # We interpret this as printing the elements of the constructed sequence A.
            print(*result_A)
