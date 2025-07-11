import sys

# Increase the recursion limit for deep calculations, although memoization helps avoid this.
sys.setrecursionlimit(2000)

# A dictionary to store results of subproblems to avoid re-calculation.
memo = {}

def get_pieces(m, k):
    """
    Recursively calculates the number of pieces for the number m / 2^k.
    m: numerator
    k: the power of 2 in the denominator (i.e., denominator is 2**k)
    """
    # If the result is already computed, return it.
    if (m, k) in memo:
        return memo[(m, k)]

    # Base case: N(0) = 0
    if m == 0:
        return 0
    
    # Base case: N(n) for integer n.
    # N(m/2^k) = N(n) if m = n * 2^k. So if k=0, it's just N(m).
    # We will only need N(1) which is 1 piece.
    if k == 0:
        return m # N(m) = m for integers

    # If the numerator m is even, simplify the fraction first.
    # N(m/2^k) = N((m/2)/(2^(k-1)))
    if m % 2 == 0:
        result = get_pieces(m // 2, k - 1)
    else:
        # If m is odd, use the recursive formula:
        # N(m/2^k) = N((m-1)/2^k) + N((m+1)/2^k) + 2
        # Note: (m-1) and (m+1) are both even.
        left_pieces = get_pieces(m - 1, k)
        right_pieces = get_pieces(m + 1, k)
        result = left_pieces + right_pieces + 2

    # Store the result before returning.
    memo[(m, k)] = result
    return result

# The number we are interested in is 13/16.
# So, m = 13 and the denominator is 16 = 2^4, so k = 4.
numerator = 13
denominator_power = 4

total_pieces = get_pieces(numerator, denominator_power)

print(f"To represent the number {numerator}/{2**denominator_power}, we can use the game {{ {numerator-1}/{2**denominator_power} | {numerator+1}/{2**denominator_power} }} = {{ {int((numerator-1)/4)}/{int(2**denominator_power/4)} | {int((numerator+1)/2)}/{int(2**denominator_power/2)} }}.")
print(f"The number of pieces for N({numerator}/{2**denominator_power}) is recursively calculated as N({(numerator-1)//4}/{2**(denominator_power-2)}) + N({(numerator+1)//2}/{2**(denominator_power-1)}) + 2.")
print(f"The total number of red or blue pieces needed is {total_pieces}.")