import numpy as np

# Step 1 & 2: Determine the quantity supplied by the producer.
# The profit function pi(Q) = 4Q + 0.03Q^3 - 0.0005Q^4 is maximized
# at the upper bound of the production constraint [0, 10].
quantity_supplied = 10

# Step 3: Determine the equilibrium price set by the producer.
# P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3
# The producer sets the price P for Q = 10.
price = 4 + 0.03 * (quantity_supplied**2) - 0.0005 * (quantity_supplied**3)

# Step 4: Calculate the total quantity demanded at the equilibrium price.
# We need to solve the cubic equation for Q_d:
# Q_d^3 - 60*Q_d^2 + 5000 = 0
# The coefficients for the polynomial Q^3 - 60Q^2 + 0Q + 5000 are [1, -60, 0, 5000].
coeffs = [1, -60, 0, 5000]
roots = np.roots(coeffs)

# The roots are the possible quantities demanded. We look for positive, real roots.
# The economically stable equilibrium is the largest positive root.
quantity_demanded = max(r.real for r in roots if r.imag == 0 and r.real > 0)

# Step 5: Calculate the excess demand.
excess_demand = quantity_demanded - quantity_supplied

print(f"The producer will supply Q_supply = {quantity_supplied} units.")
print(f"The producer will set the price at P = {price}.")
print(f"At this price, the market quantity demanded Q_d is {quantity_demanded:.4f}.")
print("\nThe final equation for excess demand is: Quantity Demanded - Quantity Supplied")
print(f"Excess Demand = {quantity_demanded:.4f} - {quantity_supplied} = {excess_demand:.4f}")
