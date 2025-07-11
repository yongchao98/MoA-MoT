import math

def combinations(n, k):
    """
    Calculates the number of combinations 'n choose k'.
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve():
    """
    Calculates the total number of different sequences (a_1,...,a_8, b_1,...,b_6)
    by computing the sum from the generating function expansion.
    """
    total_count = 0
    N_target = 48
    num_vars = 14  # 8 players in team A + 6 in team B

    print("The final answer is the sum of the following terms:")
    print("Formula: Sum of C(8,k) * C(6,j) * (-1)^(k+j) * C(61 - 7k - 9j, 13) for valid k,j")
    print("-" * 60)

    # Loop over the exponent 'k' for the (1-x^7)^8 term
    for k in range(9):  # k from 0 to 8
        # Loop over the exponent 'j' for the (1-x^9)^6 term
        for j in range(7):  # j from 0 to 6
            power = 7 * k + 9 * j
            
            if N_target - power >= 0:
                # Coefficient from (1-x^7)^8 and (1-x^9)^6 expansions
                term_coeff = combinations(8, k) * combinations(6, j) * ((-1)**(k + j))
                
                # Coefficient from (1-x)^(-14) expansion
                # C(n+d-1, d-1) for a power n
                n = N_target - power
                d = num_vars
                comb_term = combinations(n + d - 1, d - 1)
                
                term_value = term_coeff * comb_term
                total_count += term_value
                
                # Print each non-zero term in the sum
                if term_value != 0:
                    sign = "+" if term_value > 0 else "-"
                    print(f"k={k}, j={j}: {sign} {abs(term_value)}")


    print("-" * 60)
    print(f"The total number of different sequences is: {total_count}")
    
solve()