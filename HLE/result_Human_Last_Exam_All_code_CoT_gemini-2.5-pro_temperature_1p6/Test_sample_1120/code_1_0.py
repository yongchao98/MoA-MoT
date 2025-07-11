import numpy as np

# Step 1: Define the Market Demand Curve
# The individual demand is q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20.
# With 100 identical customers, Q = 100 * q_i, so q_i = Q/100.
# Substituting q_i into the equation:
# Q/100 = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
# This simplifies to the inverse demand curve P(Q):
# 100P = 400 + 3Q^2 - Q^3/20  or P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3

# Step 2: Determine the Producer's Supply Quantity (Q_s)
# Profit is Total Revenue TR(Q) = P(Q)*Q, since Marginal Cost is 0.
# TR(Q) = 4*Q + 0.03*Q^3 - 0.0005*Q^4
# Marginal Revenue MR(Q) = d(TR)/dQ = 4 + 0.09*Q^2 - 0.002*Q^3
# The producer's capacity is Q <= 10. Let's check the MR at the capacity limit.
Q_s = 10.0
marginal_revenue_at_capacity = 4 + 0.09 * (Q_s**2) - 0.002 * (Q_s**3)
print(f"The producer's marginal revenue at the capacity limit Q=10 is {marginal_revenue_at_capacity:.2f}.")
print(f"Since MR({Q_s}) = {marginal_revenue_at_capacity:.2f} is greater than the marginal cost of 0, the producer will supply at maximum capacity.")
print(f"Quantity Supplied (Q_s) = {Q_s:.2f}")
print("-" * 30)

# Step 3: Find the Equilibrium Price (P_eq)
# The producer sets the price to sell the entire supply Q_s = 10.
P_eq = 4 + 0.03 * (Q_s**2) - 0.0005 * (Q_s**3)
print(f"The equilibrium price is set by plugging Q_s = {Q_s:.2f} into the inverse demand curve P(Q).")
print(f"P_eq = 4 + 0.03*({Q_s:.2f}^2) - 0.0005*({Q_s:.2f}^3) = {P_eq:.2f}")
print("-" * 30)

# Step 4: Determine the Quantity Demanded (Q_d) at the Equilibrium Price
# At P_eq = 6.5, the demand-side equilibrium quantity Q must solve:
# 0 = 400 - 100*P_eq + 3Q^2 - Q^3/20
# Multiplying by -20 gives: Q^3 - 60Q^2 + 2000*P_eq - 8000 = 0
# Q^3 - 60Q^2 + 2000*6.5 - 8000 = 0 => Q^3 - 60Q^2 + 5000 = 0
coefficients = [1, -60, 0, 5000]
roots = np.roots(coefficients)
print(f"At P_eq = {P_eq:.2f}, we find the possible quantities demanded by solving: Q^3 - 60Q^2 + 5000 = 0.")
print(f"The possible quantities demanded are: {[f'{r:.2f}' for r in roots if np.isreal(r)]}")

# In a market with a single producer, the only rational expectation for consumers
# is the quantity the producer intends to sell, which is Q=10. This is also the
# stable equilibrium among the possible demand-side quantities.
Q_d = 10.0
print(f"The stable and rational equilibrium for quantity demanded is Q_d = {Q_d:.2f}.")
print("-" * 30)

# Step 5: Calculate Excess Demand
# Excess Demand = Quantity Demanded - Quantity Supplied
excess_demand = Q_d - Q_s
print("The final equation for excess demand is:")
print(f"Excess Demand = {Q_d:.2f} - {Q_s:.2f}")
print(f"The value of the excess demand at the equilibrium price is: {excess_demand:.2f}")