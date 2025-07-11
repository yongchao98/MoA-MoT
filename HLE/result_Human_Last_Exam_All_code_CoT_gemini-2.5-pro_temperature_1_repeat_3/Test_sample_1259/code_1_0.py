import math

def restore_sequence():
    """
    Reads a set S and attempts to restore a sequence A that could generate it.
    """
    try:
        s_input = list(map(int, input().split()))
    except (ValueError, EOFError):
        # Handle empty or invalid input
        print("-1")
        return

    if not s_input:
        print("-1")
        return

    # S is a set, so we only consider unique elements, sorted.
    s_unique_sorted = sorted(list(set(s_input)))
    
    # The smallest element of S must be the gcd of the entire sequence A.
    g = s_unique_sorted[0]
    
    # A necessary and sufficient condition is that g must divide all elements of S.
    for element in s_unique_sorted:
        if element % g != 0:
            print("-1")
            return
            
    # If the condition is met, a solution is guaranteed to exist.
    # We can construct one by interleaving the elements of S \ {g} with g.
    
    n = len(s_unique_sorted)
    
    # If S has only one element, A is just that element.
    if n == 1:
        restored_A = [g]
    else:
        # A simple construction: A = [u1, g, u2, g, ...], where ui are elements of S other than g.
        # This creates an array of length 2 * (n - 1).
        restored_A = []
        # Get elements other than g
        s_others = s_unique_sorted[1:]
        for val in s_others:
            restored_A.append(val)
            restored_A.append(g)

    # The prompt asks to "output each number in the final equation!".
    # We interpret this as printing the length of A and the elements of A.
    # However, to be cleaner, we will just print the elements of A.
    print(*restored_A)

# Execute the function
restore_sequence()