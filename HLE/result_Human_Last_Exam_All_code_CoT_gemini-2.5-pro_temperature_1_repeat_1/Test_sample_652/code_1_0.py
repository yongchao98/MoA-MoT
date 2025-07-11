import numpy as np
from scipy.optimize import fsolve

# Problem parameters
p = 0.2
wage = 20

# 1. Find the optimal q for p = 0.2
# The implicit function H(q, p) = 0 is derived in the plan.
# H(q, p) = exp(2q)(1-q) - (1.5-q) - (p/(1-p))*(1-q)^2 = 0
# We define a function for the root-finding algorithm.
def find_q_equation(q, p_val):
    return np.exp(2*q)*(1-q) - (1.5-q) - (p_val/(1-p_val))*(1-q)**2

# Solve for q using an initial guess of 0.3
q_star = fsolve(find_q_equation, 0.3, args=(p))[0]

# 2. Calculate the elasticity using the derived formula
# Elasticity_qp = (dq/dp) * (p/q)
# dq/dp = Numerator / Denominator

# Numerator of dq/dp
dq_dp_numerator = ((1 - q_star)**2) / ((1 - p)**2)

# Denominator of dq/dp
dq_dp_denominator = np.exp(2 * q_star) * (1 - 2 * q_star) + 1 + (2 * p * (1 - q_star)) / (1 - p)

# dq/dp
dq_dp = dq_dp_numerator / dq_dp_denominator

# Elasticity
elasticity = dq_dp * (p / q_star)

# 3. Output the results step-by-step
print("Step 1: Find the optimal search intensity q* for p = 0.2")
print(f"The optimal q* is: {q_star:.5f}\n")

print("Step 2: Calculate the elasticity using the formula Elasticity = (dq/dp) * (p / q*)")
print("The calculation for dq/dp = Numerator / Denominator is as follows:")
print(f"  Numerator = (1 - q*)^2 / (1 - p)^2 = {dq_dp_numerator:.5f}")
print(f"  Denominator = e^(2q*)(1-2q*) + 1 + 2p(1-q*)/(1-p) = {dq_dp_denominator:.5f}")
print(f"  dq/dp = {dq_dp_numerator:.5f} / {dq_dp_denominator:.5f} = {dq_dp:.5f}\n")

print("The final elasticity calculation is:")
print(f"  Elasticity = dq/dp * (p / q*)")
print(f"  Elasticity = {dq_dp:.5f} * ({p} / {q_star:.5f})")
print(f"  Elasticity = {elasticity:.5f}\n")

print("The elasticity of optimal q with respect to p, rounded to three decimals, is:")
print(f"{elasticity:.3f}")

# Final answer in the required format
final_answer = round(elasticity, 3)
# <<<0.221>>>