import math

# Thinking process translated into code comments and logic.

# 1. Game Setup & Interpretation
# The game is sequential: C chooses w_C, then B chooses w_B, then A chooses w_A.
# The strategic variable is `w`, but since p(w) is strictly decreasing, choosing `w`
# is equivalent to choosing a success probability `p` from (0, 1] (p=p(w)).
# Let's analyze the game in terms of the players choosing `p_C`, `p_B`, and `p_A`.
# The choices must be unique.
# Objective: Maximize the chance of winning. We model the probability of winning for agent `i` as:
# Win_i = p_i * product(1 - p_j) for j != i.
# This formula captures the idea that an agent wins if they uniquely succeed.

# 2. Backward Induction - Agent A's move
# Agent A sees p_B and p_C. A chooses p_A to maximize Win_A = p_A * (1 - p_B) * (1 - p_C).
# The term (1 - p_B) * (1 - p_C) is constant for A.
# To maximize Win_A, A must maximize p_A.
# A must choose from the available set (0, 1] excluding {p_B, p_C}.
# A's optimal strategy, p_A_star, is to choose the highest possible probability.
# - If p_B < 1 and p_C < 1, A will choose p_A = 1. This makes Win_A > 0, and Win_B = Win_C = 0.
# - If p_B = 1 or p_C = 1, then Win_A is always 0. A is indifferent to any valid p_A.

# 3. Backward Induction - Agent B's move
# Agent B sees p_C and anticipates A's rational response.
# Let's analyze B's options given C's choice, p_C.

# Case 1: C chooses p_C < 1.
# - If B chooses p_B < 1 (and p_B != p_C): A will see two probabilities less than 1. A will choose p_A = 1.
#   This results in Win_B = p_B * (1 - p_A) * (1 - p_C) = p_B * (1 - 1) * (1 - p_C) = 0.
# - If B chooses p_B = 1: A sees p_C < 1 and p_B = 1. A's winning probability is now Win_A = p_A * (1 - 1) * (1 - p_C) = 0.
#   A is indifferent. In this situation, B's winning probability is Win_B = 1 * (1 - p_A) * (1 - p_C).
#   This probability is greater than 0 as long as A doesn't choose p_A = 1 (which A cannot, as B has already chosen it).
#   Since a positive winning probability is better than a zero winning probability, B will choose p_B = 1.

# Summary of B's strategy: If p_C < 1, B's best move is to choose p_B = 1 to guarantee A cannot win and give itself a chance to win.

# 4. Backward Induction - Agent C's move
# Agent C moves first and anticipates the actions of both B and A.
# C evaluates its own options.

# Case 1: C chooses p_C < 1.
# - As determined above, B will respond by choosing p_B = 1.
# - In this scenario, C's winning probability is Win_C = p_C * (1 - p_A) * (1 - p_B) = p_C * (1 - p_A) * (1 - 1) = 0.
# - So, any choice of p_C < 1 guarantees that C will lose.

# Case 2: C chooses p_C = 1.
# - B must choose p_B < 1, since choices must be unique.
# - A must choose p_A < 1, and p_A != p_B.
# - C's winning probability is Win_C = 1 * (1 - p_A) * (1 - p_B). This value is positive.
# - B's winning probability is Win_B = p_B * (1 - p_A) * (1 - 1) = 0.
# - A's winning probability is Win_A = p_A * (1 - p_B) * (1 - 1) = 0.
#   In this situation, both B and A have a winning probability of 0 and are indifferent. But C's winning probability is positive.

# Conclusion of C's strategy: C compares getting a winning probability of 0 (by choosing p_C < 1) versus a positive winning probability (by choosing p_C = 1).
# A rational Agent C will choose p_C = 1 to maximize their chance of winning.
# This reasoning holds even under the simpler deterministic model where the agent with the highest `p` wins.

# 5. Final Calculation
# The optimal success probability for Agent C is p_C = 1.
p_C = 1.0

# The problem asks for floor(100 * p_C).
factor = 100
result_value = factor * p_C
final_answer = math.floor(result_value)

# Print the equation as requested
print(f"The optimal probability for Agent C is p_C = {p_C}.")
print("The calculation is:")
print(f"{factor} * {p_C} = {result_value}")
print(f"The floor of this value is {final_answer}.")
print("\nFinal Answer:")
print(final_answer)
