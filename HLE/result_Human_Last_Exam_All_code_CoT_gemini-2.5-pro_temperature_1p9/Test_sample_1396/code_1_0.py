import math

# This script calculates the complexity for finding a connected epsilon-envy-free
# cake allocation for 4 agents based on the most recent academic findings.

# Number of agents
n = 4

# Epsilon (ε) represents the tolerance for envy-freeness. A smaller epsilon
# means a more "fair" allocation but requires more computation (queries).
# We will use a sample value of 0.01 for this demonstration.
epsilon = 0.01

# The upper bound 'O' is given by the formula from recent research:
# O = (n^2 / ε) * log(n / ε)
# Note: 'log' refers to the natural logarithm. The base of the logarithm does not
# affect the asymptotic complexity, but we use a specific base for calculation.

# --- Step-by-Step Calculation ---

# 1. Calculate n squared
n_squared = n**2

# 2. Calculate the term inside the logarithm
log_argument = n / epsilon

# 3. Calculate the natural logarithm of the argument
log_value = math.log(log_argument)

# 4. Calculate the coefficient (n^2 / ε)
coefficient = n_squared / epsilon

# 5. Calculate the final estimated complexity value
complexity_value = coefficient * log_value


# --- Output the results with each number in the equation ---

print("The formula for the query complexity upper bound 'O' is: (n^2 / ε) * log(n / ε)")
print(f"\nFor this problem, we are given n = {n}.")
print(f"We will use a sample epsilon (ε) value of {epsilon}.")
print("-" * 60)
print("Calculating the bound step-by-step:")

# Print the equation with values substituted
print(f"\nO = ({n}^2 / {epsilon}) * log({n} / {epsilon})")

# Print the result of the first-level evaluations (n^2 and n/ε)
print(f"O = ({n_squared} / {epsilon}) * log({log_argument})")

# Print the result of the division and the logarithm calculation
print(f"O = {coefficient} * {log_value:.4f}")

# Print the final numerical result
print(f"O ≈ {complexity_value:.2f}")

print("-" * 60)
print(f"\nThus, for {n} agents and an envy tolerance of ε = {epsilon}, the number of queries")
print(f"is bounded by approximately {int(round(complexity_value))}.")
