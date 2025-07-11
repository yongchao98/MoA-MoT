import math

# Step 1: The inverse demand function is P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3.
# This is derived by aggregating individual demands:
# Q = 100 * (400 - 100P + Q/100 + 3Q^2 - Q^3/20)
# Q = 40000 - 10000P + Q + 300Q^2 - 5Q^3
# 10000P = 40000 + 300Q^2 - 5Q^3
# P = 4 + 0.03Q^2 - 0.0005Q^3

# Step 2: Determine the producer's profit-maximizing quantity.
# Profit pi(Q) = P(Q)*Q = 4Q + 0.03Q^3 - 0.0005Q^4 for Q <= 10.
# Marginal Revenue MR(Q) = d(pi)/dQ = 4 + 0.09Q^2 - 0.002Q^3.
# Marginal Cost MC(Q) = 0 for Q <= 10.
# Since MR(Q) is always positive for Q in [0, 10], the producer maximizes
# profit by supplying the maximum possible quantity.
quantity_supplied = 10

# Step 3: Find the equilibrium price P for Q_s = 10.
# P = 4 + 0.03*(10)^2 - 0.0005*(10)^3 = 4 + 3 - 0.5 = 6.5
price_equilibrium = 6.5

# Step 4: Calculate the quantity demanded at P = 6.5.
# We need to solve for Q_d in the demand equation:
# 10000 * 6.5 = 40000 + 300*(Q_d)^2 - 5*(Q_d)^3
# 65000 = 40000 + 300*(Q_d)^2 - 5*(Q_d)^3
# 5*(Q_d)^3 - 300*(Q_d)^2 + 25000 = 0
# (Q_d)^3 - 60*(Q_d)^2 + 5000 = 0
# The positive roots of this equation are Q_d = 10 and Q_d = 25 + 15*sqrt(5).
# The larger root represents a stable demand equilibrium.
quantity_demanded = 25 + 15 * math.sqrt(5)

# Step 5: Calculate the excess demand.
excess_demand = quantity_demanded - quantity_supplied

# Print the final calculation and result
print("The producer supplies Q_s = 10 units at a price P = 6.5.")
print(f"At this price, the quantity demanded is Q_d = 25 + 15 * sqrt(5) = {quantity_demanded:.4f}")
print("The excess demand is Q_d - Q_s.")
print(f"Excess Demand = (25 + 15 * sqrt(5)) - 10")
print(f"Excess Demand = 15 + 15 * sqrt(5)")
print(f"The value of the excess demand is: {excess_demand:.4f}")
