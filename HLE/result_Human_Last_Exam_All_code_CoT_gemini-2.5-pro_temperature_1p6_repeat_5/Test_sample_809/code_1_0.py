import math

def solve_and_calculate_ev():
    """
    This function performs all the calculations as described in the plan.
    1. Determines h.
    2. Solves for the probability of traversing g, P(g).
    3. Calculates the probability of traversing k, P(k).
    4. Calculates the binomial probabilities for the bet.
    5. Computes the final expected value.
    """
    # Step 1: Determine h
    # As derived, the smallest value of x for guaranteed success is 1/2.
    h = 0.5
    print(f"Step 1: The value of h is {h}\n")

    # Step 2: Solve for P(g)
    # The probabilities for the paths in g are:
    # Path 1 (direct success): h = 0.5
    # Path 3 (6-g chain): 1/6
    # Path 2 (hole): 2 * (1/6) = 1/3
    # The equation for P(g) is: P(g) = 1/2 + (1/6) * P(g)^6
    print("Step 2: Solving for P(g) using the equation P(g) = 1/2 + (1/6) * P(g)^6")
    
    # Numerical solution using fixed-point iteration
    p_g = 0.5  # Initial guess
    for i in range(10): # 10 iterations is sufficient for convergence
        p_g = 0.5 + (1/6) * (p_g**6)
    
    print(f"The calculated probability of traversing g, P(g), is: {p_g:.8f}\n")

    # Step 3: Calculate P(k)
    # P(k) is the probability of traversing 4 'g's in a chain.
    p_k = p_g**4
    print("Step 3: Calculating P(k) = P(g)^4")
    print(f"The probability of successfully traversing k, P(k), is: {p_k:.8f}\n")

    # Step 4: Calculate the probability of the opponent winning
    # This is a binomial probability calculation.
    n = 100  # number of trials
    p_success = p_k # probability of success for one trial
    
    print(f"Step 4: Calculating opponent's win probability for a binomial distribution B(n={n}, p={p_success:.8f})")
    print("Opponent wins if the number of successes (X) is less than 6 (X=0, 1, 2, 3, 4, 5).\n")

    prob_opponent_wins = 0
    print("Calculating individual probabilities P(X=i):")
    for i in range(6):
        # Calculate C(n, i) * p^i * (1-p)^(n-i)
        prob_i = math.comb(n, i) * (p_success**i) * ((1 - p_success)**(n - i))
        prob_opponent_wins += prob_i
        print(f"P(X={i}) = {prob_i:.8f}")

    print(f"\nTotal probability of opponent winning, P(X < 6), is: {prob_opponent_wins:.8f}\n")

    # Step 5: Calculate the expected value
    expected_value = 1 - 2 * prob_opponent_wins
    print("Step 5: Calculating your expected value (EV)")
    print(f"EV = 1 - 2 * P(opponent wins)")
    print(f"EV = 1 - 2 * {prob_opponent_wins:.8f}")
    print(f"EV = {expected_value:.8f}")

    # Round to the nearest cent
    final_ev = round(expected_value, 2)
    print(f"\nYour expected value of playing one game is ${final_ev:.2f}.")

# Execute the function
solve_and_calculate_ev()