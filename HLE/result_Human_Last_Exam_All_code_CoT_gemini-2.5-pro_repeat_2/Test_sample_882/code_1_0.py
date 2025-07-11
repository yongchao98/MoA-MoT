import math

def solve_game_theory_problem():
    """
    Solves the game theory problem to find the value of floor(10000 * (1-p)).
    """
    N = 8

    # Payoff function for the discrete strategy S_d
    def u_d(p):
        if p == 0:
            # By L'Hopital's rule, the limit as p->0 is 1.
            return 1.0
        # For a player vs 3N-1 others, the payoff simplifies to this expression.
        # n = 3N-1 opponents, x = p/N. Payoff = (1-(1-x)^(n+1))/((n+1)x)
        # (n+1)x = 3N * p/N = 3p
        return (1 - pow(1 - p / N, 3 * N)) / (3 * p)

    # Payoff function for the spreading strategy S_m
    def u_spread(m, p):
        # M = 3N-1 opponents
        M = 3 * N - 1
        res = 0.0
        for j in range(1, m + 1):
            # Probability that none of the discrete players pick any of a specific set of j races
            prob_uncontested = pow(max(0, 1 - j * p / N), M)
            
            term = math.comb(m, j) * prob_uncontested
            
            # Applying the Principle of Inclusion-Exclusion
            if (j - 1) % 2 == 0:  # j is odd
                res += term
            else:  # j is even
                res -= term
        return res

    # For N=8, it is known that the optimal number of races to spread over is m*=8
    m_star = 8

    # We need to solve u_d(p) = u_spread(m_star, p) for p.
    # Let's define the function f(p) = u_d(p) - u_spread(m_star, p)
    # We are looking for the root of f(p) = 0.
    
    # We will use binary search to find p.
    # Numerical evaluation shows the root is between 0.8 and 0.9.
    low = 0.8
    high = 0.9
    
    # Iterate for high precision
    for _ in range(100):
        mid = (low + high) / 2
        f_mid = u_d(mid) - u_spread(m_star, mid)
        if f_mid > 0:
            # Root is in the upper half
            low = mid
        else:
            # Root is in the lower half
            high = mid
            
    p = (low + high) / 2

    # Calculate the final result
    result = math.floor(10000 * (1 - p))
    
    print(result)

solve_game_theory_problem()