import math

def solve_game_theory_problem():
    """
    This function solves for the mixed strategy equilibrium probability 'p' in the described game
    for N=8 and then calculates the final requested value.
    """
    N = 8
    # Total players are 3N. Number of opponents for any given player is 3N-1.
    M = 3 * N - 1

    def E_discrete(p):
        """
        Calculates the expected payoff for a player using the discrete strategy.
        This is the probability of winning the single race they bet on.
        """
        if p == 0:
            # The limit of the expression as p -> 0 is 1.
            return 1.0
        
        base = 1.0 - p / N
        # The term (1-p/N)^3N
        term_pow = base**(3 * N)
        
        # The full expression for the payoff
        payoff = (1.0 - term_pow) / (3.0 * p)
        return payoff

    def E_deviation(p):
        """
        Calculates the expected payoff for a player using the uniform spread deviation.
        The payoff is the probability of winning at least one race, calculated via
        the Principle of Inclusion-Exclusion.
        """
        total_prob = 0.0
        for k in range(1, N + 1):
            # Probability of winning k specific races simultaneously. This occurs if no
            # opponent plays the discrete strategy on any of these k races.
            prob_k_wins = (1.0 - k * p / N)**M
            
            # Term from the Inclusion-Exclusion formula
            term = math.comb(N, k) * prob_k_wins
            
            if k % 2 == 1:  # Add for odd k
                total_prob += term
            else:  # Subtract for even k
                total_prob -= term
        return total_prob

    # We need to find p such that E_discrete(p) = E_deviation(p).
    # We use binary search to find the root of f(p) = E_discrete(p) - E_deviation(p).
    low = 0.0
    high = 1.0

    # Perform 100 iterations for high precision, sufficient for the required significant digits.
    for _ in range(100):
        p_guess = (low + high) / 2.0
        
        # If the search interval is no longer shrinking, break.
        if p_guess == low or p_guess == high:
            break
        
        difference = E_discrete(p_guess) - E_deviation(p_guess)
        
        # If difference > 0, p is too low, so we search in the upper half.
        # If difference < 0, p is too high, so we search in the lower half.
        if difference > 0:
            low = p_guess
        else:
            high = p_guess
            
    # The equilibrium probability p is the result of our search.
    p = high

    # The problem asks for the floor of 10000 * (1-p).
    # The numbers in this final equation are 10000, 1, and the calculated p.
    final_value = math.floor(10000.0 * (1.0 - p))
    
    print(final_value)

solve_game_theory_problem()