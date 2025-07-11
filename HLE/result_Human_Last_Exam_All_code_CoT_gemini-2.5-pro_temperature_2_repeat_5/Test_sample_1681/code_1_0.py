import numpy as np

# Payoff matrix values
R = 5  # Reward for mutual cooperation
P = 1  # Punishment for mutual defection

# The optimal quantum game sets the entanglement such that cos^2(gamma/2) = 1/3
# The formula for the new equilibrium payoff is P_eq = R*cos^2(gamma/2) + P*sin^2(gamma/2)
# which simplifies to P_eq = (R-P)*cos^2(gamma/2) + P
# With our values, this is P_eq = 4 * cos^2(gamma/2) + 1.

cos_sq_gamma_half = 1/3

# Calculate the equilibrium payoff
payoff_value = (R - P) * cos_sq_gamma_half + P

# Output the final equation with the numbers plugged in
print(f"The equilibrium payoff is calculated from the equation:")
print(f"Payoff = ({R} - {P}) * {cos_sq_gamma_half:.3f} + {P}")
print(f"Payoff = {R-P} * (1/3) + {P}")
print("Payoff = 4 * (1/3) + 1")
print("Payoff = 4/3 + 3/3")
print("Payoff = 7/3")

print("\nThe numerical value of the equilibrium payoff for each player is:")
print(payoff_value)

print("\nThe final equilibrium point is a payoff of (7/3, 7/3) for the two players.")

# We output the single value as requested by the format.
print("<<<7/3>>>")