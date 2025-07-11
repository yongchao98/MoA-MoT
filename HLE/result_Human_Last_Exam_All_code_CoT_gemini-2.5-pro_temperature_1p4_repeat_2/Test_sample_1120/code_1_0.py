import numpy as np

# Step 1 & 2: Determine the quantity supplied
# The producer's marginal cost is 0 for Q <= 10.
# The Total Revenue is TR(Q) = P(Q) * Q = (4 + 0.03*Q^2 - 0.0005*Q^3) * Q
# The Marginal Revenue is MR(Q) = 4 + 0.09*Q^2 - 0.002*Q^3.
# For Q in [0, 10], MR(Q) is always positive, so revenue is maximized at the upper bound.
q_supply = 10
print(f"The producer's profit-maximizing quantity supplied is: {q_supply}")

# Step 3: Find the equilibrium price for Q = 10
p_eq = 4 + 0.03 * (q_supply**2) - 0.0005 * (q_supply**3)
print(f"The equilibrium price set by the producer is: {p_eq}")

# Step 4 & 5: Find the total quantity demanded at the equilibrium price
# We need to solve the cubic equation: Q^3 - 60*Q^2 + 5000 = 0
# The coefficients are [1, -60, 0, 5000] for Q^3, Q^2, Q^1, Q^0
coefficients = [1, -60, 0, 5000]
roots = np.roots(coefficients)
print(f"\nThe possible quantities demanded at P={p_eq} are the roots of the equation Q^3 - 60Q^2 + 5000 = 0.")
print(f"The roots are: {roots[0]:.4f}, {roots[1]:.4f}, {roots[2]:.4f}")

# The economically meaningful positive roots are the potential demand levels.
# We choose the larger positive root as it corresponds to a stable market equilibrium.
q_demand = max(r for r in roots if r > 0)
print(f"\nThe stable market quantity demanded is the largest positive root: {q_demand:.4f}")

# Step 6: Calculate the excess demand
excess_demand = q_demand - q_supply
print(f"\nTo find the excess demand, we compute: Quantity Demanded - Quantity Supplied")
print(f"The final equation is: {q_demand:.4f} - {q_supply:.4f} = {excess_demand:.4f}")
print(f"\nThe value of the excess demand is: {excess_demand:.4f}")
