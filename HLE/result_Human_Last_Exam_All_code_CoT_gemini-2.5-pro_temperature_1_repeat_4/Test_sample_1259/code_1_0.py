import math

def solve():
    """
    This function reads a set of numbers S, determines if a valid sequence A
    can be constructed, and if so, prints one such sequence.
    For demonstration purposes, the input is hardcoded.
    """
    
    # Example input. In a real scenario, this would be read from stdin.
    # For S = {2, 6, 10}, the code would find a solution like A = [6, 2, 10, 2].
    # For S = {6, 10}, min(S)=6 but gcd(S)=2, so it would output -1.
    # For S = {6, 36, 60}, min(S)=6 and gcd(S)=6, so it finds a solution.
    S = [6, 36, 60]

    if not S:
        print(-1)
        return

    # Sort S to easily find the minimum element.
    S.sort()

    def gcd_list(numbers):
        """Computes the greatest common divisor of a list of numbers."""
        if not numbers:
            return 0  # GCD of an empty set is undefined, 0 is a common convention.
        result = numbers[0]
        for i in range(1, len(numbers)):
            result = math.gcd(result, numbers[i])
        return result

    min_s = S[0]
    gcd_s = gcd_list(S)

    # Check the necessary and sufficient condition.
    if min_s != gcd_s:
        print(-1)
    else:
        # If the condition holds, a solution is guaranteed to exist.
        # We construct one such solution.
        if len(S) == 1:
            A = S
        else:
            # A robust construction is to interleave the smallest element `g`
            # with the other elements of the set.
            # A = [s2, g, s3, g, ..., sm, g]
            A = []
            g = S[0]
            other_elements = S[1:]
            for element in other_elements:
                A.append(element)
                A.append(g)
        
        # As per the prompt, we output the "final equation", which we interpret
        # as the restored sequence A.
        # First, print the length of the sequence.
        print(len(A))
        # Then, print the elements of the sequence separated by spaces.
        print(*A)

solve()