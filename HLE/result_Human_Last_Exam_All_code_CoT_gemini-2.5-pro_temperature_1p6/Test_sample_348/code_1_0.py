import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Step 1: Define the supply and demand functions
def demand(q):
    """Demand function: P = 18 * e^(-arctan(Q))"""
    q = np.asarray(q, dtype=float)
    return 18 * np.exp(-np.arctan(q))

def supply(q):
    """Supply function: P = ln(Q^3 - 2)"""
    q = np.asarray(q, dtype=float)
    # The logarithm is only defined for Q^3 - 2 > 0.
    # We suppress errors for values outside this domain during solving.
    with np.errstate(invalid='ignore'):
        result = np.log(q**3 - 2)
    return result

# Step 2: Find the market equilibrium
def equilibrium_equation(q):
    """Equation to solve for equilibrium: Supply(Q) - Demand(Q) = 0"""
    return supply(q) - demand(q)

# An initial guess for the solver, based on preliminary analysis.
initial_guess_q = 4.5
# Use fsolve to find the root, which is the equilibrium quantity.
equilibrium_quantity_result = fsolve(equilibrium_equation, initial_guess_q)
Q_e = equilibrium_quantity_result[0]

# Calculate the equilibrium price using the demand function
P_e = demand(Q_e)

print("--- Market Equilibrium ---")
print(f"The equilibrium quantity (Q_e) is: {Q_e:.4f}")
print(f"The equilibrium price (P_e) is: {P_e:.4f}\n")


# Step 3: Calculate Consumer and Producer Surplus
# The integral of the demand function from 0 to Q_e
area_under_demand, _ = quad(demand, 0, Q_e)
consumer_surplus = area_under_demand - P_e * Q_e

# The minimum quantity for which the supply function is defined (Q^3 - 2 > 0)
q_min_supply = np.cbrt(2)
# The integral of the supply function from its starting point to Q_e
area_under_supply, _ = quad(supply, q_min_supply, Q_e)
producer_surplus = P_e * Q_e - area_under_supply

print("--- Surplus Calculation ---")
print(f"Consumer Surplus (CS) is calculated as: ∫D(Q)dQ - P_e*Q_e = {area_under_demand:.4f} - {P_e*Q_e:.4f} = {consumer_surplus:.4f}")
print(f"Producer Surplus (PS) is calculated as: P_e*Q_e - ∫S(Q)dQ = {P_e*Q_e:.4f} - {area_under_supply:.4f} = {producer_surplus:.4f}\n")


# Step 4: Calculate Total Welfare
total_welfare = consumer_surplus + producer_surplus

print("--- Total Welfare Calculation ---")
print("Total Welfare is the sum of Consumer and Producer Surplus.")
print(f"Final Equation: Total Welfare = CS + PS")
print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f}")
print(f"The final Total Welfare is: {total_welfare:.4f}")
<<<22.3934>>>