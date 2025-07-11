import math

def solve_expected_value():
    """
    This function calculates the expected value of the game based on the problem description.
    """

    # Step 1: Determine the value of h.
    # The probability of traversing f, P(f), follows P(f) = x + (1-x)P(f)^2.
    # The solutions are P(f)=1 and P(f)=x/(1-x).
    # For traversal to be guaranteed, P(f) must be 1. This requires x/(1-x) >= 1,
    # which gives x >= 0.5. The smallest such x is 0.5.
    h = 0.5

    # Step 2: Determine the success probability of g, q.
    # Path probabilities for g: P(direct)=h=0.5, P(hole)=2*P(rec), and they sum to 1.
    # 0.5 + 2*P(rec) + P(rec) = 1  =>  3*P(rec) = 0.5  =>  P(rec) = 1/6.
    # The success probability q follows the equation: q = 0.5 + (1/6) * q^6.
    # We solve this numerically through iteration.
    q = 0.5  # Initial guess
    for _ in range(20):
        q = 0.5 + (1/6) * (q**6)

    # Step 3: Determine the success probability of k, p_k.
    # k is a chain of 4 instances of g.
    p_k = q**4

    # Step 4: Calculate the opponent's win probability using the Binomial distribution.
    # n=100 trials, p=p_k success probability.
    # Opponent wins if successes X < 6 (i.e., X <= 5).
    n = 100
    p = p_k
    p_opponent_wins = 0

    # We need a combination function C(n, k) = n! / (k!(n-k)!)
    # We use memoization for factorial to be efficient.
    factorial_memo = {}
    def factorial(i):
        if i not in factorial_memo:
            factorial_memo[i] = math.factorial(i)
        return factorial_memo[i]
        
    def combinations(n_c, k_c):
        if k_c < 0 or k_c > n_c:
            return 0
        return factorial(n_c) // (factorial(k_c) * factorial(n_c - k_c))

    for i in range(6):  # Summing probabilities for i = 0, 1, 2, 3, 4, 5
        term = combinations(n, i) * (p**i) * ((1-p)**(n-i))
        p_opponent_wins += term

    # Step 5: Calculate the final expected value.
    # EV = P(You Win) * $1 + P(Opponent Wins) * (-$1)
    # EV = (1 - p_opponent_wins) - p_opponent_wins
    expected_value = 1 - 2 * p_opponent_wins
    
    print(f"The probability of one successful traversal of k is p_k = {p_k:.6f}")
    print(f"The opponent wins if the number of successes in 100 trials is less than 6.")
    print(f"The probability of the opponent winning is P(X <= 5) = {p_opponent_wins:.6f}")
    print("\nThe final expected value (EV) is calculated as:")
    print(f"EV = P(You Win) * $1 - P(Opponent Wins) * $1")
    print(f"EV = (1 - P(Opponent Wins)) - P(Opponent Wins) = 1 - 2 * P(Opponent Wins)")
    print(f"EV = 1 - 2 * {p_opponent_wins:.6f}")
    print(f"EV = {expected_value:.6f}")
    print(f"\nYour final expected value, rounded to the nearest cent, is: ${expected_value:.2f}")

solve_expected_value()