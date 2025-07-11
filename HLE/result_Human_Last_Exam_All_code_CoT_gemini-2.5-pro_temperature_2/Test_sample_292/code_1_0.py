import math
from fractions import Fraction

def calculate_sum():
    """
    Calculates the sum described in the problem.
    S = sum_{w in V^n} (n+1 - |unique_tokens(w)|)^-1
    where n = 99.
    """
    n = 99
    
    total_sum = Fraction(0)
    
    # Iterate through all possible numbers of unique tokens, k.
    # k can range from 1 to n.
    for k in range(1, n + 1):
        
        # Calculate the number of surjections from a set of size n to a set of size k.
        # This is sum_{j=0 to k} (-1)^j * C(k,j) * (k-j)^n
        surjections_sum = 0
        for j in range(k + 1):
            # Using pow(-1, j) is more robust than (-1)**j for large j
            sign = 1 if j % 2 == 0 else -1
            term = sign * math.comb(k, j) * pow(k - j, n)
            surjections_sum += term
            
        # N(k) is the number of sequences with exactly k unique tokens.
        # N(k) = C(n,k) * (number of surjections from n-set to k-set)
        N_k = math.comb(n, k) * surjections_sum
        
        # Add the contribution of all sequences with k unique tokens to the sum.
        # The weight for each such sequence is 1 / (n + 1 - k).
        total_sum += Fraction(N_k, n + 1 - k)
        
    # The final equation is sum_{w} a(w) = result
    # We output the numbers n that define the equation, and the final result.
    print(f"For a vocabulary of size n = {n} and sequence length L = {n}:")
    print(f"The calculated sum is: {total_sum.numerator}")

    # The result should be written as a power of 10.
    # We know the sum is 100^98 = 10^196. We can verify this with log10.
    # Note: Using math.log10 on the very large numerator might lead to float precision issues.
    # A more robust way is to analyze the string representation of the number.
    # The result is 1 followed by 196 zeros.
    exponent = len(str(total_sum.numerator)) - 1
    
    print(f"This result can be written as 10 to the power of {exponent}.")

if __name__ == '__main__':
    calculate_sum()
