import scipy.optimize
import scipy.stats

def solve_and_calculate_ev():
    """
    This function follows the step-by-step plan to calculate the expected value of the game.
    """
    # Step 1: Determine h
    # As derived in the plan, h is the smallest x where P(f) = 1.
    # P(f) = x / (1 - x), so P(f) >= 1 when x >= 0.5.
    h = 0.5
    print(f"Step 1: The value of h is {h}")

    # Step 2: Determine P(g), the probability of traversing g.
    # The probabilities for the paths in g are:
    # P1 (direct success) = h = 0.5
    # P3 (recursive chain) = (1 - h) / 3 = 0.5 / 3 = 1/6
    # P2 (hole) = 2 * P3 = 1/3
    # The equation for P(g) is: P_g = 0.5 + (1/6) * P_g^6
    # We need to find the root of the function f(y) = y^6 - 6y + 3 = 0.
    # The root is known to be between 0 and 1.
    g_equation = lambda y: y**6 - 6*y + 3
    p_g = scipy.optimize.brentq(g_equation, 0, 1)
    print(f"Step 2: The probability of traversing g, P(g), is {p_g:.8f}")

    # Step 3: Determine P(k), the probability of traversing k.
    # k is a chain of 4 instances of g.
    p_k = p_g ** 4
    print(f"Step 3: The probability of traversing k, P(k) = P(g)^4, is {p_k:.8f}")

    # Step 4: Calculate the expected value of the game.
    n_trials = 100
    # The opponent bets on less than 6 successful traversals (i.e., 0, 1, 2, 3, 4, or 5).
    success_threshold = 5

    # Probability of the opponent winning is P(X <= 5)
    prob_opponent_wins = scipy.stats.binom.cdf(k=success_threshold, n=n_trials, p=p_k)
    print(f"\nStep 4: Calculating Game Outcome Probabilities")
    print(f"Number of trials (n): {n_trials}")
    print(f"Probability of a single success (p): {p_k:.8f}")
    print(f"Opponent wins if successes are less than 6 (X <= {success_threshold}).")
    print(f"Probability opponent wins = P(X <= {success_threshold}) = {prob_opponent_wins:.8f}")

    # Probability of you winning is P(X >= 6) = 1 - P(X <= 5)
    prob_you_win = 1 - prob_opponent_wins
    print(f"Probability you win = P(X >= 6) = {prob_you_win:.8f}")

    # Calculate Expected Value
    # EV = (Prob you win * $1) + (Prob opponent wins * -$1)
    expected_value = (prob_you_win * 1) + (prob_opponent_wins * -1)
    
    print("\nFinal Calculation for Expected Value (EV):")
    print(f"EV = (P(win) * $1) + (P(lose) * -$1)")
    print(f"EV = ({prob_you_win:.8f} * 1) + ({prob_opponent_wins:.8f} * -1)")
    print(f"EV = {expected_value:.8f}")
    
    print(f"\nYour expected value of playing one game is: ${expected_value:.2f}")

solve_and_calculate_ev()
<<<O>>>