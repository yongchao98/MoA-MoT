import numpy as np
from scipy.optimize import root_scalar

# Step 1 & 2: Determine the producer's optimal quantity and the resulting price.
# The inverse demand curve is P(Q) = 4 + 0.03*Q**2 - 0.0005*Q**3.
# The producer's profit is pi(Q) = P(Q)*Q = 4*Q + 0.03*Q**3 - 0.0005*Q**4.
# The capacity constraint is Q <= 10.
# The marginal profit is d(pi)/dQ = 4 + 0.09*Q**2 - 0.002*Q**3.
# This marginal profit is positive for all Q in [0, 10], so profit is maximized at the constraint.
# Therefore, the quantity supplied is Q_s = 10.
Q_s = 10

# Step 3: The equilibrium price is the price that clears the market for Q_s = 10.
P_eq = 4 + 0.03 * (Q_s**2) - 0.0005 * (Q_s**3)

# Step 4: Find the demand equilibria at P_eq.
# At a given price P, the quantity demanded Q must satisfy:
# Q = 100 * (400 - 100*P + Q/100 + 3*Q**2 - Q**3/20)
# Plugging in P = 6.5, this simplifies to the cubic equation:
# Q**3 - 60*Q**2 + 5000 = 0.
# We need to find the roots of this equation.

def demand_equilibrium_eq(Q):
    """Represents the demand equilibrium equation: Q^3 - 60Q^2 + 5000 = 0."""
    return Q**3 - 60 * Q**2 + 5000

# Step 5: Identify the stable root.
# Stability analysis shows that an equilibrium Q* is stable if Q* > 40.
# The equilibrium Q=10 is unstable. The market will converge to the larger, stable root.
# We use a numerical solver to find the root in the stable region (Q > 40).
try:
    sol = root_scalar(demand_equilibrium_eq, bracket=[41, 60], method='brentq')
    Q_d = sol.root
except ValueError:
    print("Could not find a stable root in the specified bracket.")
    Q_d = float('nan')
    
# Step 6: Calculate and print the excess demand.
excess_demand = Q_d - Q_s

print("The final equation is: Excess Demand = Quantity Demanded - Quantity Supplied")
print("The numbers in this equation are:")
print(f"Quantity Demanded (stable equilibrium): {Q_d:.4f}")
print(f"Quantity Supplied (producer's choice): {Q_s:.4f}")
print(f"Excess Demand = {Q_d:.4f} - {Q_s:.4f} = {excess_demand:.4f}")
