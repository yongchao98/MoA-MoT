import math

def C(n, k):
    """
    Computes the binomial coefficient C(n, k), also known as "n choose k".
    Returns 0 if k < 0 or k > n, as is standard in combinatorics.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve():
    """
    Calculates the number of different sequences (a_1,...,a_8, b_1,...,b_6) that can occur.
    """
    total_sequences = 0
    
    # The number of sequences is the coefficient of x^48 in the expansion of
    # P(x) = (1+x+...+x^6)^8 * (1+x+...+x^8)^6.
    # This coefficient is given by the sum:
    # Sum_{i,j} (-1)^(i+j) * C(8,i) * C(6,j) * C(48-7i-9j+13, 13)
    # where the sum is over all i,j such that 7i+9j <= 48.
    
    print("The total number of sequences is the sum of the following terms:")
    print("-" * 60)

    # Iterate over possible values for i from the expansion of (1-x^7)^8
    # and j from the expansion of (1-x^9)^6
    for j in range(7):  # j can be at most 5, since 9*6 > 48
        for i in range(9): # i can be at most 6, since 7*7 > 48
            
            # This is the inclusion-exclusion constraint
            if 7 * i + 9 * j > 48:
                continue
            
            # Binomial coefficients from (1-x^7)^8 and (1-x^9)^6
            c1 = C(8, i)
            c2 = C(6, j)
            
            # The exponent for the (1-x)^-14 term
            k = 48 - 7 * i - 9 * j
            c3 = C(k + 13, 13)
            
            sign = (-1)**(i + j)
            
            term = sign * c1 * c2 * c3
            
            # We only print non-zero terms to keep the output clean.
            if term != 0:
                op = "+" if sign > 0 else "-"
                print(f"Term (i={i}, j={j}): {op} C(8,{i})*C(6,{j})*C({k+13},13) = {term}")
                
            total_sequences += term
            
    print("-" * 60)
    print(f"The final sum, which is the total number of possible sequences, is:")
    print(total_sequences)
    print(f"\n<<< {total_sequences} >>>")

solve()