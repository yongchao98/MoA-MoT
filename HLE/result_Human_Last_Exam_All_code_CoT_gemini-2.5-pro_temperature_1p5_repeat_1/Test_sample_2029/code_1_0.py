import math

def solve_match_combinations():
    """
    This function calculates the number of different score sequences possible in the tournament.
    The method uses generating functions to model the problem, leading to a summation formula
    for the coefficient of x^48 in a polynomial expansion.
    """

    # The problem is equivalent to finding the coefficient of x^48 in the expansion of:
    # P(x) = (1 + x + ... + x^6)^8 * (1 + x + ... + x^8)^6
    #
    # This can be rewritten using geometric series formula as:
    # P(x) = ((1-x^7)/(1-x))^8 * ((1-x^9)/(1-x))^6
    # P(x) = (1-x^7)^8 * (1-x^9)^6 * (1-x)^(-14)
    #
    # The final count is the coefficient of x^48, which we find by summing up terms.
    # The general formula derived is:
    # Sum_{i=0 to 8} Sum_{j=0 to 6} [ C(8, i) * C(6, j) * (-1)^(i+j) * C(48 - (7*i + 9*j) + 14 - 1, 14 - 1) ]
    #
    # Let's explain the numbers in this final formula:
    print("The number of different score sequences is calculated based on the following equation derived from generating functions.")
    print("The numbers in the equation represent the parameters of the problem:")
    print("   8: Number of players in team A.")
    print("   6: Number of players in team B.")
    print("   7: Maximum wins for a player in team A + 1 (from the term 1-x^7).")
    print("   9: Maximum wins for a player in team B + 1 (from the term 1-x^9).")
    print("  48: Total number of games played.")
    print("  14: Total number of players (8 + 6), which becomes the exponent in (1-x)^(-14).")
    print("  13: Is 14 - 1, used in the combination formula C(n+k-1, k-1).")
    print("  61: Is 48 + 13, which simplifies the formula C(48 - k + 13, 13) to C(61-k, 13).")
    
    print("\nThe formula is:")
    print("Total = Sum_{i=0 to 8} Sum_{j=0 to 6} [ (-1)^(i+j) * C(8, i) * C(6, j) * C(61 - 7*i - 9*j, 13) ]\n")

    total_count = 0
    
    # Loop for the expansion of (1-x^7)^8, where i is the power of the (-x^7) term
    for i in range(8 + 1):
        # Loop for the expansion of (1-x^9)^6, where j is the power of the (-x^9) term
        for j in range(6 + 1):
            
            # The exponent of x contributed by the numerator terms is (7*i + 9*j)
            # The final coefficient comes from the term x^(48 - (7*i + 9*j)) in the expansion of (1-x)^(-14)
            # This is given by C( (48 - (7*i + 9*j)) + 14 - 1, 14 - 1 )
            
            n = 61 - (7 * i + 9 * j)
            k = 13
            
            # math.comb(n, k) is 0 if n < k, so we only add non-zero terms
            if n >= k:
                sign = (-1)**(i + j)
                
                # C(8, i) from the first binomial expansion
                term1 = math.comb(8, i)
                
                # C(6, j) from the second binomial expansion
                term2 = math.comb(6, j)
                
                # C(n, k) from the (1-x)^(-14) expansion
                term3 = math.comb(n, k)
                
                total_count += sign * term1 * term2 * term3
                
    print(f"The total number of different sequences is: {total_count}")

solve_match_combinations()