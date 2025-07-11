import numpy as np

# Step 1 & 2: Determine producer's supply
# The producer's profit is TR = P*Q. The inverse demand is P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3.
# Total Revenue TR(Q) = 4Q + 0.03*Q^3 - 0.0005*Q^4.
# Marginal Revenue MR(Q) = d(TR)/dQ = 4 + 0.09*Q^2 - 0.002*Q^3.
# For Q in [0, 10], MR(Q) is always positive (e.g., MR(0)=4, MR(10)=11).
# Since MC=0 and MR>0, the producer maximizes profit by producing at maximum capacity.
quantity_supplied = 10.0

# Step 3: Find the equilibrium price P* for Q_s = 10
# From the aggregate demand equation: 100P = 400 + 3Q^2 - Q^3/20
price_star = (400 + 3 * quantity_supplied**2 - (quantity_supplied**3) / 20) / 100

print(f"The producer will supply Q_s = {quantity_supplied} units.")
print(f"The equilibrium price P* set by the producer is ${price_star:.2f}.")
print("-" * 30)

# Step 4: Find the quantity demanded (Q_d) at P*
# We need to solve the cubic equation for Q_d: 0 = 400 - 100*P* + 3*Q_d^2 - Q_d^3/20
# This simplifies to: 0.05*Q_d^3 - 3*Q_d^2 + (100*P* - 400) = 0
# 0.05*Q_d^3 - 3*Q_d^2 - 250 = 0 -> Correction: 0 = -250 + 3Q_d^2 - 0.05Q_d^3 => 0.05Q_d^3 - 3Q_d^2 + 250 = 0
coefficients = [0.05, -3, 0, 250]
roots = np.roots(coefficients)

# Step 5: Filter for positive real solutions and choose the stable equilibrium
q_demand_solutions = [r.real for r in roots if abs(r.imag) < 1e-9 and r.real > 0]

# The demand function has two equilibria at this price. One is Q=10, which is an unstable
# equilibrium. The other higher value is the stable equilibrium, which represents the
# actual market demand.
quantity_demanded = max(q_demand_solutions)

print(f"At the price of ${price_star:.2f}, there are two possible demand levels: {q_demand_solutions[0]:.2f} and {q_demand_solutions[1]:.2f}.")
print(f"The higher value, Q_d = {quantity_demanded:.2f}, represents the stable demand equilibrium.")
print("-" * 30)


# Step 6: Calculate excess demand
excess_demand = quantity_demanded - quantity_supplied

print("Excess demand is the difference between the quantity demanded and the quantity supplied at the equilibrium price.")
print("Final Equation:")
print(f"{quantity_demanded:.2f} (Quantity Demanded) - {quantity_supplied:.2f} (Quantity Supplied) = {excess_demand:.2f}")
