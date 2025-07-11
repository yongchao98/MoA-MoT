import numpy as np
import math

# Step 1: Define problem parameters
Q_supply = 10
N_customers = 100

# Step 2: Find the equilibrium price (P_eq) for the producer to sell Q_supply = 10 units.
# The market demand equation is derived from the individual demand function:
# q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
# and the market equilibrium condition Q = 100 * q_i.
# This gives the relation: 0 = 400 - 100P + 3Q^2 - Q^3/20.
# To sell Q = 10, the producer must set a price P such that this equation holds.
# 100P = 400 + 3*(10**2) - (10**3)/20
# 100P = 400 + 3*100 - 1000/20
# 100P = 400 + 300 - 50
# 100P = 650
P_eq = 650 / 100

print(f"The producer's profit-maximizing quantity to sell is its capacity, Q_supply = {Q_supply}.")
print(f"To sell this quantity, the producer sets a price P where the market demand equals 10.")
print(f"Solving the demand equation for P when Q=10 gives the equilibrium price: P_eq = {P_eq}")

# Step 3: Find all possible quantities demanded at the equilibrium price P_eq.
# Substitute P_eq = 6.5 back into the demand relation:
# 0 = 400 - 100*(6.5) + 3Q^2 - Q^3/20
# 0 = 400 - 650 + 3Q^2 - Q^3/20
# 0 = -250 + 3Q^2 - Q^3/20
# Multiply by -20 to get a standard polynomial form:
# 0 = 5000 - 60Q^2 + Q^3
# We solve the cubic equation: Q^3 - 60Q^2 + 0Q + 5000 = 0

# Coefficients for the polynomial solver
coeffs = [1, -60, 0, 5000]
roots = np.roots(coeffs)

# Step 4: Identify the relevant demand quantities.
# We are interested in positive real roots as they represent valid market quantities.
positive_real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-9 and r.real > 0])

# The roots are the possible quantities demanded. One is Q=10. The other is the high-demand equilibrium.
Q_d_low = positive_real_roots[0]
Q_d_high = positive_real_roots[1]

print(f"\nAt the price P = {P_eq}, there are multiple possible levels of total market demand.")
print("This is due to the network effects in the demand function.")
print(f"Solving the cubic equation Q^3 - 60Q^2 + 5000 = 0 gives the possible demand levels:")
print(f"Possible Demand 1: {Q_d_low:.4f}")
print(f"Possible Demand 2: {Q_d_high:.4f}")

# Step 5: Calculate the excess demand.
# The question asks for the value of excess demand. This refers to the situation where
# the market coordinates on the high-demand equilibrium.
# Excess Demand = Higher Quantity Demanded - Quantity Supplied
excess_demand = Q_d_high - Q_supply

print("\nThe 'excess demand' refers to the difference between the high-demand equilibrium and the quantity supplied.")
print("\nThe final equation is: Excess Demand = High Quantity Demanded - Quantity Supplied")
print(f"Excess Demand = {Q_d_high} - {Q_supply}")
print(f"\nThe calculated value of the excess demand is: {excess_demand}")