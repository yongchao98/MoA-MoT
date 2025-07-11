import math
try:
    from scipy.stats import binom
except ImportError:
    print("SciPy library not found. Please install it using: pip install scipy")
    exit()

def solve_task():
    """
    Solves the multi-step probability problem and calculates the expected value of the game.
    """

    # --- Step 1: Find h ---
    # The probability of traversing f, P(f), is given by the recursive formula:
    # P(f) = x * 1 + (1 - x) * P(f)^2
    # Rearranging this gives a quadratic equation for P(f): (1-x)P(f)^2 - P(f) + x = 0.
    # The solutions are P(f) = 1 and P(f) = x/(1-x).
    # For success to be "guaranteed", P(f) must be 1. This becomes the only stable
    # solution when the other solution is no longer a valid probability less than 1,
    # i.e., x/(1-x) >= 1, which implies 2x >= 1, so x >= 0.5.
    # The smallest value of x for this is 0.5.
    h = 0.5

    print("Step 1: The probability of the direct path in f, x, must be at least 0.5 for the traversal to be guaranteed.")
    print(f"h = {h}")
    print("-" * 30)

    # --- Step 2: Find P(g) ---
    # In g, we have three paths with probabilities p1, p2, p3.
    # p1 = h = 0.5
    # p2 = 2 * p3
    # p1 + p2 + p3 = 1 => 0.5 + 2*p3 + p3 = 1 => 3*p3 = 0.5 => p3 = 1/6.
    # p2 = 2 * (1/6) = 1/3.
    p_chain_g = (1 - h) / 3.0

    # The probability of traversing g, P(g), is given by:
    # P(g) = p1*1 + p2*0 + p3*P(g)^6  => P(g) = h + p_chain_g * P(g)^6
    # We solve this equation numerically.
    pg = 0.5  # Initial guess
    for _ in range(10): # Iterate to find the fixed point
        pg = h + p_chain_g * (pg ** 6)

    print("Step 2: The probability of successfully traversing system g, P(g), is calculated.")
    print(f"The recursive equation is: P(g) = {h} + {p_chain_g:.4f} * P(g)^6")
    print(f"Solving numerically gives P(g) = {pg:.8f}")
    print("-" * 30)

    # --- Step 3: Find P(k) ---
    # k is a chain of 4 instances of g.
    num_g_in_k = 4
    pk = pg ** num_g_in_k
    
    print("Step 3: The probability of successfully traversing system k, P(k), is calculated.")
    print(f"k is a chain of {num_g_in_k} g's, so P(k) = P(g)^{num_g_in_k}.")
    print(f"P(k) = ({pg:.8f})^{num_g_in_k} = {pk:.8f}")
    print("-" * 30)

    # --- Step 4: Calculate the Expected Value ---
    # The game is 100 independent traversals of k. This is a binomial distribution.
    # n = 100 trials, p = P(k) is the success probability.
    # Opponent wins if number of successes X is < 6 (i.e., X <= 5).
    # We need P(X <= 5).
    n_trials = 100
    win_cutoff = 5  # Opponent wins if successes are 0, 1, 2, 3, 4, or 5.
    
    p_opponent_wins = binom.cdf(win_cutoff, n_trials, pk)
    
    # Expected Value = (+1) * P(You Win) + (-1) * P(Opponent Wins)
    # EV = P(X >= 6) - P(X <= 5)
    # EV = (1 - P(X <= 5)) - P(X <= 5)
    # EV = 1 - 2 * P(X <= 5)
    
    val_1 = 1
    val_2 = 2
    ev = val_1 - val_2 * p_opponent_wins
    
    print("Step 4: The expected value (EV) of the game is calculated.")
    print(f"The probability of the opponent winning is P(Successes <= {win_cutoff}) for a Binomial({n_trials}, {pk:.4f}) distribution.")
    print(f"The final EV equation is: EV = {val_1} - {val_2} * P(Successes <= {win_cutoff})")
    print(f"Each number in the final equation: {val_1}, {val_2}, {p_opponent_wins:.8f}")
    print("-" * 30)
    
    print(f"Final calculation: EV = {val_1} - {val_2} * {p_opponent_wins:.8f} = {ev:.8f}")
    print(f"Your expected value of playing one game is ${ev:.2f}.")

solve_task()