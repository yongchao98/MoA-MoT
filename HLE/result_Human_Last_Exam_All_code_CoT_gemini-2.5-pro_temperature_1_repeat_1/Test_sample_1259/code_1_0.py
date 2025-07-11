import math
from functools import reduce

def solve():
    """
    This function restores a sequence A from its given set of subarray GCDs, S.

    Plan:
    1. First, check if a solution can exist. A solution exists if and only if
       the minimum element of S is equal to the greatest common divisor (GCD) of all elements in S.
    2. If this condition fails, report that no solution exists.
    3. If the condition holds, construct a valid sequence A. A provably correct method is to
       create a sequence by interleaving the minimum element 'g' with the other elements of S.
       For instance, if S = {g, s1, s2}, a valid A is [g, s1, g, s2].
    4. Print the restored sequence A.
    """

    # For demonstration, a sample input S is used.
    # Replace this with actual input as needed.
    S = [2, 4, 6, 12]
    
    print(f"Given the set S = {S}")

    # Helper function to compute GCD of a list of numbers
    def gcd_list(numbers):
        if not numbers:
            return 0
        return reduce(math.gcd, numbers)

    if not S:
        print("Result: -1")
        print("Reason: The input set S cannot be empty.")
        return

    # Use a set to handle potential duplicates and sort for consistent order
    s_list = sorted(list(set(S)))
    m = len(s_list)

    # Check the necessary and sufficient condition
    g = s_list[0]
    overall_gcd = gcd_list(s_list)

    if g != overall_gcd:
        print("Result: -1")
        print(f"Reason: No valid construction exists. The minimum element ({g}) is not the GCD of the set ({overall_gcd}).")
        return

    # If the condition holds, construct a valid sequence A
    if m == 1:
        A = s_list
    else:
        # Construction: A = [g, s1, g, s2, ...]
        A = []
        s_other = s_list[1:]
        for s_val in s_other:
            A.append(g)
            A.append(s_val)

    print(f"A valid sequence A can be constructed with length {len(A)}.")
    # The prompt asks to "output each number in the final equation!"
    # We interpret this as printing the restored sequence A.
    print("Restored sequence A:", end=" ")
    # The following line prints the numbers of sequence A separated by spaces
    print(*A)

# Execute the solution
solve()