import math
import numpy as np

# Step 1: Define the relationship between price P and total quantity Q
# Individual demand: q_i = 400 - 100*P + Q/100 + 3*Q^2 - Q^3/20
# Total demand: Q = 100 * q_i
# Q = 100 * (400 - 100*P + Q/100 + 3*Q^2 - Q^3/20)
# Q/100 = 400 - 100*P + Q/100 + 3*Q^2 - Q^3/20
# 0 = 400 - 100*P + 3*Q^2 - Q^3/20
# 100*P = 400 + 3*Q^2 - Q^3/20
# P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3
# This is the inverse demand curve.

# The producer maximizes profit Pi(Q) = P(Q)*Q for Q in [0, 10].
# Since MC=0, Profit = Revenue = 4*Q + 0.03*Q^3 - 0.0005*Q^4.
# The derivative dPi/dQ = 4 + 0.09*Q^2 - 0.002*Q^3 is always positive on [0, 10].
# So, the profit is maximized at the capacity limit.
Q_supply = 10.0

# Step 2: Calculate the equilibrium price set by the producer for Q_supply = 10
P_eq = 4 + 0.03 * (Q_supply**2) - 0.0005 * (Q_supply**3)

# Step 3: Find the quantity demanded (Q_D) at this price P_eq.
# We must solve the demand equation for Q_D at P = P_eq = 6.5.
# Q_D = 100 * (400 - 100*6.5 + Q_D/100 + 3*Q_D^2 - Q_D^3/20)
# This simplifies to the cubic equation: Q_D^3 - 60*Q_D^2 + 5000 = 0.
cubic_coeffs = [1, -60, 0, 5000]
roots = np.roots(cubic_coeffs)

# Step 4: Identify the stable demand equilibrium.
# Stability analysis shows an equilibrium Q* is stable if Q* > 40.
stable_Q_D = 0
for r in roots:
    # We only consider positive real roots
    if r.imag == 0 and r.real > 0:
        if r.real > 40:
            stable_Q_D = r.real
            break

# Step 5: Calculate and print the excess demand.
excess_demand = stable_Q_D - Q_supply

# For printing the final equation, we use the analytical solution for the roots
# The stable root is 25 + 15 * sqrt(5)
# The quantity supplied is 10
q_d_val = 25 + 15 * math.sqrt(5)
q_s_val = 10
final_val = 15 + 15 * math.sqrt(5)

print("The quantity supplied by the producer is Q_S = 10.")
print(f"The equilibrium price set by the producer is P = {P_eq}.")
print("The quantity demanded by consumers at this price, Q_D, is the stable root of the equation Q^3 - 60*Q^2 + 5000 = 0.")
print(f"The stable quantity demanded is Q_D = 25 + 15 * sqrt(5) ≈ {q_d_val:.4f}.")
print("\nExcess demand is the quantity demanded minus the quantity supplied.")
print("The final equation for excess demand is:")
print(f"(25 + 15 * sqrt(5)) - {q_s_val} = 15 + 15 * sqrt(5)")
print(f"\nThe value of the excess demand is approximately {final_val:.4f}.")