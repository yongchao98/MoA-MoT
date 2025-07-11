import numpy as np

# Step 1 & 2: Determine producer's optimal supply
# The producer's profit function is increasing for Q in [0, 10].
# Therefore, the producer supplies the maximum possible quantity.
Q_supply = 10.0

# Step 3: Calculate the equilibrium price
# The price is determined by substituting Q_supply into the inverse demand curve:
# P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3
P_eq = 4 + 0.03 * (Q_supply**2) - 0.0005 * (Q_supply**3)

print(f"The producer's optimal quantity to supply is {Q_supply:.2f} units.")
print(f"The resulting equilibrium price set by the producer is ${P_eq:.2f}.")

# Step 4: Set up the equation for quantity demanded at P_eq
# We need to solve the cubic equation: Q^3 - 60Q^2 + 5000 = 0
# The coefficients are [1, -60, 0, 5000] for Q^3, Q^2, Q^1, and the constant term.
coefficients = [1, -60, 0, 5000]

# Step 5: Solve for Q_demand and calculate excess demand
# Find the roots of the cubic equation
roots = np.roots(coefficients)

# The quantity demanded is the largest positive real root, which represents the
# stable market equilibrium demand.
Q_demand = max(r.real for r in roots if np.isreal(r) and r.real > 0)

print(f"At the price of ${P_eq:.2f}, the total quantity demanded by the market is {Q_demand:.2f} units.")

# Calculate the excess demand
excess_demand = Q_demand - Q_supply

print("\nThe final calculation for excess demand is:")
print(f"Excess Demand = Quantity Demanded - Quantity Supplied")
print(f"Excess Demand = {Q_demand:.2f} - {Q_supply:.2f} = {excess_demand:.2f}")
