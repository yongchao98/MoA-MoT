import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Step 1: Define the supply and demand functions
def supply(Q_S):
    """Supply function P = ln(Q_S^3 - 2)"""
    # The function is defined only for Q_S^3 - 2 > 0
    with np.errstate(invalid='ignore'): # Suppress warnings for values outside domain
        if np.any(Q_S**3 - 2 <= 0):
            return np.inf # Return infinity if outside the domain for the solver
        return np.log(Q_S**3 - 2)

def demand(Q_D):
    """Demand function P = 18 * e^(-arctan(Q_D))"""
    return 18 * np.exp(-np.arctan(Q_D))

# Step 2: Find the market equilibrium
def equilibrium_equation(Q):
    """Equation to solve for equilibrium: Demand(Q) - Supply(Q) = 0"""
    return demand(Q) - supply(Q)

# Initial guess for the solver. Must be in the domain of the supply function, i.e., Q > 2^(1/3) ≈ 1.26.
initial_guess = 2.0
# Use fsolve to find the equilibrium quantity (Q_e)
Q_e, = fsolve(equilibrium_equation, initial_guess)
# Calculate the equilibrium price (P_e) using the demand function
P_e = demand(Q_e)

print(f"--- Market Equilibrium ---")
print(f"Equilibrium Quantity (Qe): {Q_e:.4f}")
print(f"Equilibrium Price (Pe): {P_e:.4f}\n")

# Step 3: Calculate Consumer Surplus (CS)
# CS = ∫[0 to Qe] Demand(Q) dQ - (Pe * Qe)
integral_demand, _ = quad(demand, 0, Q_e)
consumer_surplus = integral_demand - (P_e * Q_e)

# Step 4: Calculate Producer Surplus (PS)
# The supply curve starts where Q^3 - 2 > 0, so Q > 2^(1/3)
q_min_supply = (2)**(1/3)
# PS = (Pe * Qe) - ∫[q_min_supply to Qe] Supply(Q) dQ
integral_supply, _ = quad(supply, q_min_supply, Q_e)
producer_surplus = (P_e * Q_e) - integral_supply

# Step 5: Calculate Total Welfare (TW)
total_welfare = consumer_surplus + producer_surplus

print(f"--- Welfare Calculation ---")
print(f"Consumer Surplus = ∫D(Q)dQ - Pe*Qe = {integral_demand:.4f} - {P_e * Q_e:.4f} = {consumer_surplus:.4f}")
print(f"Producer Surplus = Pe*Qe - ∫S(Q)dQ = {P_e * Q_e:.4f} - {integral_supply:.4f} = {producer_surplus:.4f}\n")
print(f"--- Final Result ---")
print(f"Total Welfare = Consumer Surplus + Producer Surplus")
print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f} = {total_welfare:.4f}")

<<<22.5029>>>