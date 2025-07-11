import numpy as np
import math

# Step 1 & 2: Determine the producer's optimal quantity to supply.
# The inverse demand curve is P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3.
# The producer's profit function is pi(Q) = P(Q)*Q for Q <= 10.
# The profit is pi(Q) = 4*Q + 0.03*Q^3 - 0.0005*Q^4.
# The marginal revenue is MR(Q) = 4 + 0.09*Q^2 - 0.002*Q^3.
# Since the marginal cost is 0 and MR(Q) is positive for all Q in [0, 10],
# the producer's profit is maximized by supplying the maximum possible quantity.
Qs = 10.0
print(f"The profit-maximizing quantity supplied by the producer is Q_s = {Qs}")

# Step 3: Find the equilibrium price set by the producer.
# The price is set to sell all Q_s = 10 units.
# P_eq = P(10) = 4 + 0.03*(10)^2 - 0.0005*(10)^3
P_eq = 4 + 0.03 * (Qs**2) - 0.0005 * (Qs**3)
print(f"The equilibrium price set by the producer is P_eq = {P_eq}")

# Step 4: Calculate the total quantity demanded (Qd) at the equilibrium price.
# At P_eq = 6.5, the quantity demanded Qd is a solution to the cubic equation:
# Qd^3 - 60*Qd^2 + 5000 = 0
coeffs = [1, -60, 0, 5000]
roots = np.roots(coeffs)

# We identify the positive roots, as quantity must be positive.
positive_real_roots = [r.real for r in roots if np.isreal(r) and r.real >= 0]

# One root is the market-clearing quantity, 10. The other positive root
# represents another possible equilibrium for demand.
# To find the non-zero excess demand, we select the larger root.
Qd = max(positive_real_roots)
print(f"At the equilibrium price, the quantity demanded is Qd = {Qd:.4f}")

# Step 5: Calculate and print the excess demand.
excess_demand = Qd - Qs

print("\nThe final equation for excess demand is: Excess Demand = Qd - Qs")
print(f"Plugging in the numbers: Excess Demand = {Qd:.4f} - {Qs:.1f}")
print(f"The value of the excess demand is: {excess_demand:.4f}")

# For verification, the exact analytical solution for the larger quantity demanded is 25 + 15*sqrt(5),
# leading to an excess demand of 15 + 15*sqrt(5).
exact_excess_demand = 15 + 15 * math.sqrt(5)
# print(f"Exact value for verification: {exact_excess_demand:.4f}")

final_answer = excess_demand
print(f"\nFinal Answer: {final_answer}")
<<<48.5410>>>