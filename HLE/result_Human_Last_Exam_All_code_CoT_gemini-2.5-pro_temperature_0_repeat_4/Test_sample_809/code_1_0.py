import math

def solve_g_prob():
    """
    Solves the equation p_g^6 - 6*p_g + 3 = 0 for p_g using Newton's method.
    This equation defines the probability of successfully traversing structure g.
    """
    # Start with an initial guess. From analysis, the root is close to 0.5.
    p_g = 0.5
    # Iterate a few times to converge to the root. 10 iterations is sufficient.
    for _ in range(10):
        f_p_g = p_g**6 - 6 * p_g + 3
        f_prime_p_g = 6 * p_g**5 - 6
        # Avoid division by zero, though unlikely here as p_g=1 is not the root.
        if f_prime_p_g == 0:
            break
        p_g = p_g - f_p_g / f_prime_p_g
    return p_g

# Step 1: Determine h
# The probability of traversing f, P(f), follows P(f) = x + (1-x)P(f)^2.
# This traversal is guaranteed (P(f)=1) for x >= 1/2.
# The smallest such value is h.
h = 0.5

# Step 2: Determine P(g), the probability of traversing g.
# The probabilities for the paths in g are P(direct)=h=0.5, P(hole)=1/3, P(chain)=1/6.
# The probability of traversing g, p_g, follows: p_g = 0.5 * 1 + (1/6) * p_g^6
# We solve the equation p_g^6 - 6*p_g + 3 = 0 for p_g.
p_g = solve_g_prob()

# Step 3: Determine P(k), the probability of traversing k.
# k is a chain of 4 g's, so P(k) = P(g)^4.
p_k = p_g**4

# Step 4: Calculate the probability of the opponent winning.
# This is a binomial probability problem with n=100 trials and p=p_k success probability.
# The opponent wins if the number of successes, S, is less than 6 (i.e., S <= 5).
n = 100
p = p_k
prob_opponent_wins = 0
for i in range(6):  # Sum probabilities for S = 0, 1, 2, 3, 4, 5
    term = math.comb(n, i) * (p**i) * ((1 - p)**(n - i))
    prob_opponent_wins += term

# Step 5: Calculate your expected value.
# EV = (+1) * P(you win) + (-1) * P(you lose)
# EV = (1 - prob_opponent_wins) - prob_opponent_wins = 1 - 2 * prob_opponent_wins
expected_value = 1 - 2 * prob_opponent_wins

# --- Output the results ---
# The final equation for the expected value is EV = 1 - 2 * P(S < 6).
# Here are the numbers that go into that calculation:
print(f"The value of h is: {h}")
print(f"The probability of traversing g, P(g), is: {p_g}")
print(f"The probability of traversing k, P(k), is: {p}")
print(f"The number of trials is: {n}")
print(f"The probability of your opponent winning (S < 6) is: {prob_opponent_wins}")
print(f"The final equation for your expected value is: $1.00 - 2 * {prob_opponent_wins:.6f}")
print(f"Your expected value is: ${expected_value:.6f}")
print(f"Rounded to the nearest cent, your expected value is: ${expected_value:.2f}")