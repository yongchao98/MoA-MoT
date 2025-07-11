import math

def solve_probability():
    """
    This function solves the entire problem step-by-step as outlined in the plan.
    """
    # Step 1: Find h
    # The probability of traversing f, P(f), satisfies P(f) = x + (1-x) * P(f)^2.
    # The smallest value of x for which P(f)=1 is guaranteed is the point where the
    # stable fixed point of the recurrence becomes 1. This occurs at x = 0.5.
    h = 0.5

    # Step 2: Define parameters for g and find P(g)
    # g has three paths with probabilities p1, p2, p3.
    # p1 = h = 0.5
    # p1 + p2 + p3 = 1 and p2 = 2 * p3
    # 0.5 + 2*p3 + p3 = 1  =>  3*p3 = 0.5  => p3 = 1/6
    # p2 = 2 * (1/6) = 1/3
    p_g_direct = h
    p_g_chain = 1.0 / 6.0
    
    # The success probability P(g) satisfies the equation:
    # P(g) = p_g_direct * 1 (success) + p_g_chain * P(g)^6
    # Let y = P(g). We need to solve y = 0.5 + (1/6)y^6, which is y^6 - 6y + 3 = 0.
    # We solve this numerically using a bisection search.
    def g_equation(y):
        return y**6 - 6*y + 3

    low, high = 0.0, 1.0
    # 100 iterations are sufficient for high precision.
    for _ in range(100):
        mid = (low + high) / 2
        if g_equation(mid) * g_equation(low) > 0:
            low = mid
        else:
            high = mid
    p_g = (low + high) / 2

    # Step 3: Find P(k), the success probability of system k
    # k is a chain of 4 instances of g.
    num_g_in_k = 4
    p_k = p_g ** num_g_in_k

    # Step 4: Calculate the probability of the opponent winning
    # This is a binomial distribution problem: n=100 trials, p=p_k probability of success.
    # The opponent wins if the number of successes X is < 6 (i.e., X <= 5).
    n = 100
    k_max = 5
    p_opponent_wins = 0
    for i in range(k_max + 1):
        # Binomial probability: C(n, i) * p^i * (1-p)^(n-i)
        term = math.comb(n, i) * (p_k**i) * ((1-p_k)**(n-i))
        p_opponent_wins += term

    # Step 5: Calculate the expected value of the game
    # EV = ($1 * P(You Win)) + (-$1 * P(Opponent Wins))
    # P(You Win) = 1 - p_opponent_wins
    # EV = (1 - p_opponent_wins) - p_opponent_wins = 1 - 2 * p_opponent_wins
    expected_value = 1 - 2 * p_opponent_wins

    # Print the final equations with the calculated numbers
    print("--- Final Equation Calculation ---")
    print(f"1. The probability of successfully traversing a single system k is p = P(g)^{num_g_in_k}.")
    print(f"   P(g) is the solution to P(g) = {p_g_direct} + {p_g_chain:.4f} * P(g)^6, which gives P(g) = {p_g:.6f}.")
    print(f"   Therefore, p = ({p_g:.6f})^{num_g_in_k} = {p_k:.6f}.")
    print()
    print(f"2. The probability of you losing (opponent winning) is P(lose) = P(successes < 6).")
    print(f"   For n={n} trials and success probability p={p_k:.6f}, this is the sum of binomial probabilities:")
    print(f"   P(lose) = Sum from i=0 to {k_max} of [ C({n}, i) * ({p_k:.6f})^i * (1-{p_k:.6f})^({n}-i) ]")
    print(f"   The result is P(lose) = {p_opponent_wins:.6f}.")
    print()
    print(f"3. The final expected value is EV = $1 * (1 - P(lose)) - $1 * P(lose) = 1 - 2 * P(lose).")
    print(f"   EV = 1 - 2 * {p_opponent_wins:.6f} = {expected_value:.6f}")
    print("---")
    print(f"Your expected value of playing one game, rounded to the nearest cent, is: ${expected_value:.2f}")


solve_probability()