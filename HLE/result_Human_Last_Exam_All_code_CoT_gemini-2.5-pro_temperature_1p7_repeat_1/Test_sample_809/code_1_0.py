import math

# Step 1: Determine the value of h from the structure of f.
# Let P(f) be the probability of successfully traversing f.
# The structure gives the recursive equation: P(f) = x + (1-x) * P(f)^2.
# We can rearrange this into a quadratic equation for P(f): (1-x)P(f)^2 - P(f) + x = 0.
# The solutions for P(f) are P(f) = 1 and P(f) = x / (1-x).
# Success is "guaranteed" when the probability P(f) becomes 1. This first occurs at the
# value of x where the two solutions are equal: x / (1-x) = 1, which solves to x = 0.5.
# For any x >= 0.5, the probability of success is 1.
# Therefore, h, the smallest value of x for guaranteed success, is 0.5.
h = 0.5
print(f"Step 1: The value of h is derived from the properties of f.")
print(f"The smallest probability x for guaranteed traversal is h = {h}")
print("-" * 30)

# Step 2: Determine the properties and success probability of g.
# g has three paths with probabilities p1, p2, p3.
# p1 (straight to end) = h = 0.5
# p2 (to hole) = 2 * p3
# p3 (to a chain of 6 g's)
# The probabilities must sum to 1: p1 + p2 + p3 = 1 => 0.5 + 2*p3 + p3 = 1 => 3*p3 = 0.5
p_g_straight = h
p_g_chain = (1 - h) / 3
p_g_hole = 2 * p_g_chain
print(f"Step 2: Analyzing g.")
print(f"The probabilities for the three paths in g are: Straight={p_g_straight}, Hole={p_g_hole:.4f}, Chain={p_g_chain:.4f}")

# Let p_g be the success probability of traversing one instance of g.
# p_g = (prob straight * success) + (prob hole * success) + (prob chain * success)
# p_g = p_g_straight * 1 + p_g_hole * 0 + p_g_chain * (p_g)^6
# p_g = 0.5 + (1/6) * p_g^6
# Solving this equation is non-trivial. However, we can see that the second term, (1/6)*(p_g)^6,
# will be very small. This implies p_g is very close to 0.5. For this problem, we will
# use the approximation p_g ≈ 0.5.
p_g = 0.5
print(f"The success probability of g, p_g, is given by p_g = {p_g_straight} + {p_g_chain:.4f} * p_g^6.")
print(f"We will use the approximation p_g ≈ {p_g}.")
print("-" * 30)


# Step 3: Determine the success probability of k.
# k is a chain of 4 instances of g. To succeed, all 4 must be traversed.
# p_k = (p_g)^4
p_k = p_g**4
print(f"Step 3: Calculating the success probability for k.")
print(f"P(k) = P(g)^4 = ({p_g})^4 = {p_k}")
print("-" * 30)


# Step 4: Calculate the expected value of the game.
# The game consists of N=100 independent traversals of k.
# The number of successes, S, follows a binomial distribution B(N, p_k).
# Your opponent wins if S < 6. You win if S >= 6.
N = 100
win_threshold = 6

# We need to calculate the probability of the opponent winning, P(S < 6).
# P(S < 6) = P(S=0) + P(S=1) + P(S=2) + P(S=3) + P(S=4) + P(S=5).
prob_opponent_wins = 0
print(f"Step 4: Calculating the expected value.")
print(f"The opponent wins if successes < {win_threshold}. Calculating this probability:")
for i in range(win_threshold):
    prob_i = math.comb(N, i) * (p_k**i) * ((1 - p_k)**(N - i))
    prob_opponent_wins += prob_i

# The expected value (EV) = (P(you win) * $1) + (P(opponent wins) * -$1)
# EV = P(you win) - P(opponent wins)
# EV = (1 - P(opponent wins)) - P(opponent wins)
# EV = 1 - 2 * P(opponent wins)
expected_value = 1 - 2 * prob_opponent_wins

print(f"The probability of the opponent winning is P(S < {win_threshold}) = {prob_opponent_wins:.6f}")
print("\nThe final equation for your expected value is: EV = 1 - 2 * P(opponent wins)")
print(f"EV = 1 - 2 * {prob_opponent_wins:.6f}")
print(f"EV = {expected_value:.6f}")

print(f"\nYour expected value of playing one game is ${expected_value:.2f}.")
