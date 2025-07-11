import math

def solve_for_p_g():
    """
    Numerically solves the equation p_g^6 - 6*p_g + 3 = 0 for p_g in [0, 1]
    using the bisection method.
    """
    def g_equation(p):
        return p**6 - 6 * p + 3

    low = 0.5
    high = 1.0
    # Iterate for a fixed number of steps for sufficient precision
    for _ in range(100):
        mid = (low + high) / 2
        if g_equation(mid) > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2

def main():
    """
    Main function to execute all steps and calculate the expected value.
    """
    # Step 1: Define h
    h = 0.5

    # Step 2: Solve for P(g), the probability of traversing g
    # P(g) = 0.5 + (1/6) * P(g)^6  => P(g)^6 - 6*P(g) + 3 = 0
    p_g = solve_for_p_g()

    # Step 3: Find P(k), the probability of traversing k
    # k is a chain of 4 instances of g
    p_k = p_g**4

    # Step 4: Calculate the probability of the opponent winning the bet
    # Opponent wins if the number of successes X in n=100 trials is less than 6.
    n = 100
    threshold = 6
    
    # Calculate P(X < 6) which is Sum(P(X=i) for i from 0 to 5)
    # P(X=i) is the binomial probability mass function
    prob_opponent_wins = 0.0
    for i in range(threshold):
        try:
            # P(X=i) = C(n, i) * p^i * (1-p)^(n-i)
            term = math.comb(n, i) * (p_k**i) * ((1 - p_k)**(n - i))
            prob_opponent_wins += term
        except (ValueError, OverflowError):
            # This handles potential precision issues, though unlikely with these numbers
            continue

    # Step 5: Calculate the expected value
    prob_you_win = 1.0 - prob_opponent_wins
    
    # E = ($1 * P(You Win)) - ($1 * P(You Lose))
    expected_value = 1.0 * prob_you_win - 1.0 * prob_opponent_wins

    # Print the numbers used in the final equation and the final result
    print(f"The probability of successfully traversing a single instance of g, P(g), is: {p_g:.6f}")
    print(f"The probability of successfully traversing a single instance of k, P(k), is: {p_k:.6f}")
    print("-" * 20)
    print("Expected Value Calculation:")
    print(f"The probability of you winning (>= 6 successes) is: {prob_you_win:.4f}")
    print(f"The probability of your opponent winning (< 6 successes) is: {prob_opponent_wins:.4f}")
    print(f"Expected Value = ($1.00 * {prob_you_win:.4f}) - ($1.00 * {prob_opponent_wins:.4f})")
    print(f"Your expected value of playing one game is: ${expected_value:.2f}")

if __name__ == "__main__":
    main()