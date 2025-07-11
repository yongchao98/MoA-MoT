import math
import sys

def solve():
    """
    Reads a set S from stdin, and if possible, constructs a sequence A
    such that the set of GCDs of all contiguous subarrays of A is S.
    """
    try:
        # Read the number of elements in S from the first line.
        # This line is often present in competitive programming inputs.
        # The problem description on some platforms provides it. Let's handle it.
        # If no count is provided, we can deduce it from the next line.
        line = sys.stdin.readline()
        if not line.strip(): # Handle empty input
             return
        # Assuming the first line might just be the count, we primarily
        # care about the list of numbers on the next line.
        # A more robust parser would handle both cases. For now, let's
        # assume one line of space-separated numbers.
        try:
             s_list_input = list(map(int, line.strip().split()))
             # In case the first number was 'm', we slice it off.
             if len(s_list_input) > 1 and s_list_input[0] == len(s_list_input) - 1:
                 s_list_input = s_list_input[1:]

        except (IOError, ValueError):
             print("-1")
             return

    except (IOError, ValueError):
        print("-1") # Or handle as appropriate
        return
    
    # The input is a set S. Let's form the canonical sorted list from it.
    S_set = set(s_list_input)
    s_list = sorted(list(S_set))

    if not s_list:
        print(-1)
        return

    # Define a helper for gcd of a list of numbers.
    def gcd_of_list(numbers):
        if not numbers:
            return 0
        result = numbers[0]
        for i in range(1, len(numbers)):
            result = math.gcd(result, numbers[i])
        return result

    # Check the necessary and sufficient condition:
    # min(S) must be equal to gcd(all elements of S).
    g = s_list[0]
    gcd_s = gcd_of_list(s_list)
    
    if g != gcd_s:
        print(-1)
    else:
        # If the condition holds, a solution is guaranteed to exist.
        # Construct a sequence A by interleaving g with the other elements of S.
        
        # Handle the trivial case where S has only one element.
        if len(s_list) == 1:
            A = [g]
        else:
            # Let s_prime be the list of elements in S other than g.
            s_prime = s_list[1:]
            
            # Construct A = [g, s'_1, g, s'_2, ...]. This always works.
            A = []
            for element in s_prime:
                A.append(g)
                A.append(element)
        
        # Print the length of the constructed sequence.
        print(len(A))
        
        # Print the elements of the sequence, space-separated.
        # The problem mentions "final equation", we interpret this as printing the sequence.
        print(*A)

solve()