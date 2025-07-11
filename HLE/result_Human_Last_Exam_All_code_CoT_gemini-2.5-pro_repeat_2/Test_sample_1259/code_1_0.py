import math
from functools import reduce

def restore_A_from_S(S_input: list[int]):
    """
    Given a set S, attempts to restore a sequence A that could generate it.
    Prints the length and elements of A if successful, otherwise prints -1.
    """
    if not S_input:
        print("-1")
        return

    # Create a sorted list of unique elements from the input set S
    s_list = sorted(list(set(S_input)))
    n_s = len(s_list)

    # Helper function to compute GCD of a list of numbers
    def gcd_list(numbers):
        if not numbers:
            return 0
        return reduce(math.gcd, numbers)

    # 1. Check the necessary and sufficient condition: min(S) == gcd(S)
    min_val = s_list[0]
    gcd_of_s = gcd_list(s_list)

    if min_val != gcd_of_s:
        print("-1")
        return

    # 2. If the condition holds, a solution is guaranteed. Construct one.
    # A simple, guaranteed construction is to interleave the non-minimum elements
    # with the minimum element 'g'.
    g = min_val
    
    # If S has only one element {g}, then A is just [g].
    if n_s == 1:
        A = [g]
    else:
        A = []
        # Get all elements other than the minimum one.
        s_others = s_list[1:]
        for s_other in s_others:
            A.append(s_other)
            A.append(g)
        # For S={s1,s2,s3} with s1=g, this produces A=[s2,g,s3,g].
        # This is a valid construction, though not always the shortest.
    
    # Print the "final equation" as the length of A and the sequence A itself.
    print(len(A))
    print(*A)

# Example usage with the counterexample for option I:
# S = {6, 12, 18}, where a solution shorter than |S| exists.
# Our generic constructor will find a valid, but not the shortest, A.
print("Example for S = {6, 12, 18}:")
restore_A_from_S([6, 12, 18])

# Example usage where the shortest A is longer than S
print("\nExample for S = {6, 60, 84, 210}:")
restore_A_from_S([6, 60, 84, 210])

# Example usage where no solution exists
print("\nExample for S = {2, 3, 6}:")
restore_A_from_S([2, 3, 6])