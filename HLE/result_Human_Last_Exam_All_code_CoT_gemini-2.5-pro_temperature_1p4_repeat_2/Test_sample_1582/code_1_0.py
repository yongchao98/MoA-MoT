import sys

# This script demonstrates that a Markov chain satisfying the given conditions
# is not positive recurrent by calculating a lower bound for its escape probability.

# --- 1. Define an example Markov Chain and the function f ---

# The Markov chain is a random walk on Z with a drift away from the origin.
# For x > 0, p(x, x+1) = 0.6, p(x, x-1) = 0.4
# For x < 0, p(x, x-1) = 0.6, p(x, x+1) = 0.4
# This chain is irreducible.

# Let the finite set be A = {-1, 0, 1}.
A = {-1, 0, 1}

# Let f(x) = |x|. This function is non-negative and f(x) -> infinity as x -> infinity.
def f(x):
  return abs(x)

# We check the submartingale condition: E[f(X_1) | X_0=x] - f(x) >= 0 for x not in A.
# For x > 1 (e.g., x is not in A):
# E[f(X_1)|X_0=x] - f(x) = (0.6 * f(x+1) + 0.4 * f(x-1)) - f(x)
#                       = (0.6 * (x+1) + 0.4 * (x-1)) - x
#                       = (0.6x + 0.6 + 0.4x - 0.4) - x
#                       = (x + 0.2) - x = 0.2 >= 0.
# The condition holds. A symmetric calculation shows it holds for x < -1.

# --- 2. Set up parameters for the escape probability calculation ---

# According to the proof, the probability of escaping to infinity (P_esc) is bounded by:
# P_esc >= (f(x0) - M) / (C - M)

# Let's choose the parameters.

# M must be greater than any f(a) for a in A.
# f(A) = {f(-1), f(0), f(1)} = {1, 0, 1}.
# max(f(A)) is 1. Let's choose M.
M = 10.0
print(f"Step 1: Define a region B_M = {{x | f(x) <= M}}.")
print(f"We choose M = {M}, which is greater than max(f(a)) for a in A = {max(f(x) for x in A)}.")
print("-" * 20)

# x0 is a starting point outside B_M, i.e., f(x0) > M.
x0 = 15
f_x0 = f(x0)
print(f"Step 2: Choose a starting state x0 such that f(x0) > M.")
print(f"We choose x0 = {x0}. f(x0) = {f_x0}, which is greater than M = {M}.")
print("-" * 20)

# C is an arbitrary constant larger than f(x0).
C = 100.0
print(f"Step 3: Choose a large constant C > f(x0).")
print(f"We choose C = {C}, which is greater than f(x0) = {f_x0}.")
print("-" * 20)

# --- 3. Calculate and print the lower bound for the escape probability ---

# The bound is (f(x0) - M) / (C - M)
numerator = f_x0 - M
denominator = C - M
escape_prob_bound = numerator / denominator

print(f"Step 4: Calculate the lower bound for the escape probability.")
print(f"The formula is: P(escape) >= (f(x0) - M) / (C - M)")
# Using sys.stdout.write to prevent the final print from adding a newline
# before the final answer tag.
sys.stdout.write(f"Plugging in the numbers:\n")
sys.stdout.write(f"P(escape) >= ({f_x0} - {M}) / ({C} - {M})\n")
sys.stdout.write(f"P(escape) >= {numerator} / {denominator}\n")
sys.stdout.write(f"P(escape) >= {escape_prob_bound}\n")
print("-" * 20)
print(f"Since the probability of escape is bounded above zero ({escape_prob_bound:.4f} > 0),")
print("the chain is transient and therefore cannot be positive recurrent.")
