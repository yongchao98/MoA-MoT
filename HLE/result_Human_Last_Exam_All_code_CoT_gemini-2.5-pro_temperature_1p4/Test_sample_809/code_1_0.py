import math

def solve_for_p_g():
    """
    Solves the equation y = 1/2 + (1/6) * y^6 for P(g) using fixed-point iteration.
    """
    # Start with an initial guess. Since p1=1/2, the value should be slightly above 1/2.
    y = 0.5
    # Iterate to find the fixed point. 20 iterations is sufficient for high precision.
    for _ in range(20):
        y = 0.5 + (1/6) * (y**6)
    return y

def binomial_probability(n, k, p):
    """
    Calculates the probability of k successes in n trials.
    P(X=k) = C(n, k) * p^k * (1-p)^(n-k)
    """
    if p < 0 or p > 1:
        return 0
    # math.comb(n, k) calculates "n choose k"
    return math.comb(n, k) * (p**k) * ((1-p)**(n-k))

def main():
    """
    Main function to calculate the expected value of the game.
    """
    # Step 1: h is 1/2. p1 for g is h.
    h = 0.5
    
    # Step 2: Solve for P(g)
    p_g = solve_for_p_g()

    # Step 3: Calculate the success probability of a single trial of k
    # k is a chain of 4 instances of g
    p_k_trial = p_g**4
    
    # Step 4: Define the parameters for the binomial distribution
    n_trials = 100
    # The bet is that the number of successes is less than 6
    k_threshold = 6
    
    # Step 5: Calculate the probability of the opponent winning
    # Opponent wins if successes are 0, 1, 2, 3, 4, or 5.
    p_opponent_wins = 0
    
    print("Calculating the opponent's win probability P(X < 6):")
    # This loop calculates P(X=0) + P(X=1) + ... + P(X=5)
    for i in range(k_threshold):
        prob_i = binomial_probability(n_trials, i, p_k_trial)
        print(f"  P(X={i}) = {prob_i:.6f}")
        p_opponent_wins += prob_i
        
    p_you_win = 1 - p_opponent_wins

    # Step 6: Calculate the expected value
    # E = (1 * P(You win)) + (-1 * P(Opponent wins))
    expected_value = p_you_win - p_opponent_wins

    print("\n--- Final Calculation ---")
    print(f"Probability h = {h}")
    print(f"Probability of traversing g, P(g) = {p_g:.6f}")
    print(f"Probability of traversing k, P(k) = P(g)^4 = {p_k_trial:.6f}")
    print(f"Number of trials = {n_trials}")
    print(f"Opponent wins if successes < {k_threshold}")
    
    print("\n--- Equation for Expected Value ---")
    # Output the numbers in the final equation as requested
    win_payout = 1
    loss_payout = -1
    print(f"P(You Win) = {p_you_win:.6f}")
    print(f"P(Opponent Wins) = {p_opponent_wins:.6f}")
    print(f"Expected Value = (P(You Win) * ${win_payout}) + (P(Opponent Wins) * ${loss_payout})")
    print(f"Expected Value = ({p_you_win:.6f} * ${win_payout}) + ({p_opponent_wins:.6f} * ${loss_payout}) = ${expected_value:.4f}")
    
    print(f"\nYour expected value of playing one game is ${expected_value:.2f}.")

if __name__ == "__main__":
    main()
<<<O>>>