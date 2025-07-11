import numpy as np
from scipy.optimize import brentq

# Step 1 & 2: Define the Marginal Revenue (MR) function.
# The monopolist maximizes profit where MR = MC. Since MC = 0, we need to find the root of MR(Q) = 0.
# TR(Q) = P(Q) * Q = (4 + 0.03*Q^2 - 0.005*Q^3) * Q = 4Q + 0.03*Q^3 - 0.005*Q^4
# MR(Q) = d(TR)/dQ = 4 + 0.09*Q^2 - 0.02*Q^3
def marginal_revenue(q):
    """Calculates the marginal revenue for a given quantity q."""
    return 4 + 0.09 * q**2 - 0.02 * q**3

# Step 3: Solve for the equilibrium quantity.
# The producer's capacity is 10. We need to find the profit-maximizing quantity Q* by solving MR(Q)=0.
# We can see that MR(0) = 4 and MR(10) = 4 + 0.09*100 - 0.02*1000 = 4 + 9 - 20 = -7.
# Since the function is continuous and changes sign between 0 and 10, a root must exist in this interval.
try:
    # Use a numerical solver (Brent's method) to find the root in the interval (0, 10).
    q_star = brentq(marginal_revenue, 0, 10)
except (RuntimeError, ValueError) as e:
    print(f"Could not find root: {e}")
    q_star = None

if q_star is not None:
    # The optimal quantity Q* is within the capacity limit (Q* < 10).
    # Therefore, the producer will supply this quantity.
    quantity_supplied = q_star
    
    # At the profit-maximizing price P(Q*), the quantity demanded by the market is, by definition, Q*.
    quantity_demanded = q_star
    
    # Step 4: Calculate and print the excess demand.
    excess_demand = quantity_demanded - quantity_supplied
    
    # The final equation is Excess Demand = Quantity Demanded - Quantity Supplied
    print("The final equation for excess demand is:")
    print(f"Excess Demand = {quantity_demanded:.10f} - {quantity_supplied:.10f}")
    print("\nResult:")
    print(f"The value of the excess demand is {excess_demand:.10f}")
