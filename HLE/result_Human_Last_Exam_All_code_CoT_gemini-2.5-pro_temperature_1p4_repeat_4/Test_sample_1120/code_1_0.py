import numpy as np
import math

# Step 1 & 2: Producer's Profit Maximization
# The market demand equation is derived from q_i = Q/100:
# Q/100 = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
# 0 = 400 - 100P + 3Q^2 - Q^3/20
# 100P = 400 + 3Q^2 - Q^3/20
# P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3
#
# The producer's Total Revenue is TR(Q) = P(Q) * Q for Q in [0, 10]
# TR(Q) = 4Q + 0.03*Q^3 - 0.0005*Q^4
# To find the maximum, we check the marginal revenue MR(Q) = d(TR)/dQ
# MR(Q) = 4 + 0.09*Q^2 - 0.002*Q^3
# For Q in [0, 10], MR(Q) is always positive (MR(0)=4, MR(10)=11).
# Since MR > 0, TR(Q) is always increasing on [0, 10].
# Therefore, the producer maximizes revenue by selling the maximum possible quantity.
q_supply = 10.0
print(f"The producer's profit-maximizing quantity to supply is Q_supply = {q_supply}")

# Step 3: Find the Equilibrium Price
# The producer sets the price to sell exactly 10 units.
# P_eq = 4 + 0.03*(10)^2 - 0.0005*(10)^3
p_eq = 4 + 0.03 * (q_supply**2) - 0.0005 * (q_supply**3)
print(f"The producer will set the equilibrium price to P_eq = {p_eq}")

# Step 4: Calculate Total Quantity Demanded
# At P = 6.5, we find Q_demand by solving:
# 6.5 = 4 + 0.03*Q^2 - 0.0005*Q^3
# This simplifies to the cubic equation: Q^3 - 60Q^2 + 5000 = 0
coefficients = [1, -60, 0, 5000]
roots = np.roots(coefficients)

# Filter for positive, real roots
positive_real_roots = [r.real for r in roots if r.imag == 0 and r.real > 0]
print(f"At P = {p_eq}, the possible demand quantities are: {', '.join([f'{r:.4f}' for r in positive_real_roots])}")

# Step 5: Identify the Stable Demand Equilibrium
# The demand curve is stable where its slope is negative.
# Slope dP/dQ = 0.06*Q - 0.0015*Q^2
def check_stability(q_val):
    slope = 0.06 * q_val - 0.0015 * (q_val**2)
    return slope < 0

# We assume the market settles at the stable equilibrium
q_demand = 0
for r in positive_real_roots:
    if check_stability(r):
        q_demand = r
        print(f"The stable demand equilibrium is Q_demand = {q_demand:.4f}")

# Step 6: Calculate Excess Demand
excess_demand = q_demand - q_supply

print("\n--- Final Calculation ---")
print("The final equation for excess demand is: Q_demand - Q_supply")
print(f"{q_demand:.4f} - {q_supply:.4f} = {excess_demand:.4f}")

# Final Answer
# The exact value is 15 + 15 * sqrt(5)
exact_value = 15 + 15 * math.sqrt(5)
print(f"\nThe exact value of the excess demand is 15 + 15*sqrt(5), which is approximately {exact_value:.4f}.")