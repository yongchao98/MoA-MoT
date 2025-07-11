import math

def solve():
    """
    This function counts the number of "allowed" pairs (a, b) with 1 <= a, b <= 1000.
    The problem of determining if a pair is "allowed" simplifies to checking if the
    greatest common divisor (GCD) of a and b is 1.
    """

    # The upper limit for a and b.
    limit = 1000
    
    # Initialize a counter for the number of allowed pairs.
    allowed_pairs_count = 0
    
    # Iterate through all possible pairs (a, b) within the specified range.
    for a in range(1, limit + 1):
        for b in range(1, limit + 1):
            # The pair (a,b) is allowed if and only if gcd(a,b) = 1.
            if math.gcd(a, b) == 1:
                allowed_pairs_count += 1
    
    # Print the total count of allowed pairs.
    print(allowed_pairs_count)

solve()