import numpy as np
from scipy.optimize import brentq

# 1. Define the inverse demand curve P(Q) based on the problem statement.
# P(Q) = 4 + 0.03*Q^2 - 0.005*Q^3
def p_q(q):
    """Calculates price for a given quantity Q."""
    return 4 + 0.03 * q**2 - 0.005 * q**3

# 2. Define the Marginal Revenue curve MR(Q).
# MR(Q) = d(P(Q)*Q)/dQ = 4 + 0.09*Q^2 - 0.02*Q^3
def mr_q(q):
    """Calculates marginal revenue for a given quantity Q."""
    return 4 + 0.09 * q**2 - 0.02 * q**3

# 3. Find the profit-maximizing quantity (Q*) where MR(Q) = MC.
# Since MC = 0, we solve for MR(Q) = 0.
# The producer's capacity is 10. We solve for Q in the interval [0, 10].
try:
    # Use a numerical solver (Brent's method) to find the root of MR(Q) = 0.
    # We found that the root is between Q=7 and Q=8, so we search in that interval.
    q_equilibrium = brentq(mr_q, 7, 8)
except ValueError:
    # Handle case where no root is found in the interval, though we know one exists.
    q_equilibrium = None

# 4. Calculate excess demand at the equilibrium.
# At the equilibrium set by the monopolist, quantity supplied (Qs) equals
# the profit-maximizing quantity.
qs = q_equilibrium

# The equilibrium price is P(q_equilibrium). At this price, the quantity
# demanded (Qd) is, by definition of the demand curve, q_equilibrium.
qd = q_equilibrium

# Excess demand is Qd - Qs.
excess_demand = qd - qs

# 5. Print the final equation.
if q_equilibrium is not None:
    print(f"The producer's profit-maximizing quantity is {q_equilibrium:.4f} units.")
    p_equilibrium = p_q(q_equilibrium)
    print(f"The producer sets the equilibrium price to {p_equilibrium:.4f}.")
    print("\nAt this equilibrium:")
    print(f"Quantity Demanded (Qd) = {qd:.4f}")
    print(f"Quantity Supplied (Qs) = {qs:.4f}")
    print("\nThe excess demand is calculated as Qd - Qs:")
    print(f"{qd:.4f} - {qs:.4f} = {excess_demand:.4f}")
else:
    print("Could not determine the equilibrium quantity.")
