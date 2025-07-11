import math

# Given parameters
alpha = 2.5
gamma = 2.1
epsilon = 0.05
confidence = 0.99
delta = 1 - confidence

# The design effect captures the sampling complexity due to graph structure.
# A smaller gamma implies a more skewed degree distribution (harder sampling).
# A larger alpha for similarity implies a less heavy tail (easier sampling).
# We model the design effect as a ratio of the exponents' differences from 1.
deff_num = (gamma - 1)
deff_den = (alpha - 1)
deff = deff_num / deff_den

# The statistical requirement is derived from a Hoeffding-like bound.
# We combine it with the design effect to find the required ratio r.
# The formula is r = Deff * ln(2/δ) * ε²
ln_term = math.log(2 / delta)
eps_squared_term = epsilon**2

r = deff * ln_term * eps_squared_term

# --- Output the results ---
print("Calculation Steps:")
print(f"1. Confidence level = {confidence}, so δ = {delta:.2f}")
print(f"2. Tolerance ε = {epsilon}")
print(f"3. Power-law exponent γ = {gamma}")
print(f"4. Pareto shape parameter α = {alpha}")
print("")
print("The formula for the sampling ratio 'r' is derived by combining a design effect with statistical requirements:")
print("r = [(γ - 1) / (α - 1)] * ln(2 / δ) * ε²")
print("")
print("Plugging in the values:")
print(f"r = [({gamma} - 1) / ({alpha} - 1)] * ln(2 / {delta:.2f}) * {epsilon}²")
print(f"r = [{deff_num} / {deff_den}] * ln({2/delta}) * {eps_squared_term}")
print(f"r = {deff:.4f} * {ln_term:.4f} * {eps_squared_term}")
print(f"r = {r:.4f}")
