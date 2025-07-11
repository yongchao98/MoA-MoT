import numpy as np
import math

# Step 1: Determine the producer's profit-maximizing quantity supplied (Q_supply).
# The producer has a capacity of 10 units with a marginal cost of 0.
# We first need the inverse demand curve P(Q) that the producer faces.
# The total market demand Q is the sum of 100 identical individual demands q_i.
# Q = 100 * q_i, so q_i = Q / 100.
# The individual demand is: q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
# Substituting q_i = Q/100, we establish the market-level relationship between P and Q:
# Q/100 = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
# 0 = 400 - 100P + 3Q^2 - Q^3/20
# From this, we derive the inverse demand function P(Q):
# 100P = 400 + 3Q^2 - Q^3/20  =>  P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3

# The producer's Total Revenue is TR(Q) = P(Q) * Q = 4Q + 0.03*Q^3 - 0.0005*Q^4.
# The Marginal Revenue is MR(Q) = d(TR)/dQ = 4 + 0.09*Q^2 - 0.002*Q^3.
# The producer's Marginal Cost is MC = 0 for Q <= 10.
# At Q=0, MR(0) = 4. At Q=10, MR(10) = 4 + 0.09*(100) - 0.002*(1000) = 4 + 9 - 2 = 11.
# Since MR is positive for all Q between 0 and 10, the producer's profit increases
# with each unit sold up to the capacity limit.
# Thus, the producer supplies their maximum possible quantity.
Q_supply = 10

# Step 2: Calculate the equilibrium price (P_eq) set by the producer.
# The price is determined by the demand for the 10 units supplied.
P_eq = 4 + 0.03 * (Q_supply**2) - 0.0005 * (Q_supply**3)

# Step 3: Determine the total quantity demanded (Q_demand) at the equilibrium price.
# We plug P_eq back into the market demand relationship:
# 0 = 40000 - 10000*P + 300*Q^2 - 5*Q^3
# 0 = 40000 - 10000*(6.5) + 300*Q^2 - 5*Q^3
# 0 = -25000 + 300*Q^2 - 5*Q^3
# This simplifies to the cubic equation: Q^3 - 60*Q^2 + 5000 = 0.

# We solve this cubic equation to find the possible values for Q_demand.
coeffs = [1, -60, 0, 5000]
roots = np.roots(coeffs)

# The roots are the possible quantities demanded. We need the stable, positive, real root.
# The demand curve is upward sloping at Q=10 (unstable) and downward sloping for the larger root (stable).
# We therefore select the largest positive root as the actual quantity demanded.
Q_demand = 0
for r in roots:
    if np.isreal(r) and r > Q_demand:
        Q_demand = np.real(r)

# Step 4: Calculate the excess demand.
excess_demand = Q_demand - Q_supply

# Print the final result and the components of the final equation.
print(f"The quantity supplied by the producer is Q_supply = {Q_supply}.")
print(f"The price set by the producer is P_eq = {P_eq}.")
print(f"At this price, the market's stable quantity demanded is Q_demand = {Q_demand:.4f}.")
print("\nThe excess demand is calculated as Q_demand - Q_supply.")
print(f"Final Equation: {Q_demand:.4f} - {Q_supply} = {excess_demand:.4f}")
print(f"\nThe value of the excess demand is {excess_demand:.4f}.")

# The exact symbolic answer for the excess demand is 15 + 15*sqrt(5)
exact_answer = 15 + 15 * math.sqrt(5)
print(f"The exact value is 15 + 15 * sqrt(5) ≈ {exact_answer:.10f}")
print(f"<<<{exact_answer}>>>")