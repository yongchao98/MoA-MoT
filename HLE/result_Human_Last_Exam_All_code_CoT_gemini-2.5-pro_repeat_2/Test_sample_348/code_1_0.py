import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Step 1: Define the supply and demand functions
def demand(Q):
    """Demand function: P = 18 * e^(-arctan(Q))"""
    return 18 * np.exp(-np.arctan(Q))

def supply(Q):
    """Supply function: P = ln(Q^3 - 2)"""
    # The function is defined for Q^3 - 2 > 0.
    # We handle the domain for the numerical integrator.
    return np.log(Q**3 - 2)

# Step 2: Find the market equilibrium
def equilibrium_equation(Q):
    """Equation to find equilibrium, where Demand(Q) - Supply(Q) = 0."""
    # We need to handle the domain of the supply function for the solver.
    # If Q is outside the domain, return a large value to push the solver away.
    if Q**3 - 2 <= 0:
        return 1e10
    return demand(Q) - supply(Q)

# The supply function is defined for Q > 2^(1/3) ≈ 1.26.
# We'll use an initial guess greater than this, e.g., 4.0.
initial_guess_Q = 4.0
Q_e = fsolve(equilibrium_equation, initial_guess_Q)[0]
P_e = demand(Q_e)

print("--- Market Equilibrium ---")
print(f"The equilibrium quantity (Q_e) is: {Q_e:.4f}")
print(f"The equilibrium price (P_e) is: {P_e:.4f}\n")

# Step 3: Calculate Consumer Surplus (CS)
# CS = ∫[0 to Q_e] D(Q) dQ - (P_e * Q_e)
integral_demand, _ = quad(demand, 0, Q_e)
consumer_surplus = integral_demand - (P_e * Q_e)

print("--- Consumer Surplus (CS) Calculation ---")
print(f"Integral of Demand from 0 to {Q_e:.4f} = {integral_demand:.4f}")
print(f"Market Revenue (P_e * Q_e) = {P_e * Q_e:.4f}")
print(f"Consumer Surplus = {integral_demand:.4f} - {P_e * Q_e:.4f} = {consumer_surplus:.4f}\n")

# Step 4: Calculate Producer Surplus (PS)
# The supply curve is only defined for Q > 2^(1/3). This is our lower integration bound.
Q_min_supply = 2**(1/3)
# PS = (P_e * Q_e) - ∫[Q_min to Q_e] S(Q) dQ
integral_supply, _ = quad(supply, Q_min_supply, Q_e)
producer_surplus = (P_e * Q_e) - integral_supply

print("--- Producer Surplus (PS) Calculation ---")
print(f"Supply curve starts at Q_min = {Q_min_supply:.4f}")
print(f"Market Revenue (P_e * Q_e) = {P_e * Q_e:.4f}")
print(f"Integral of Supply from {Q_min_supply:.4f} to {Q_e:.4f} = {integral_supply:.4f}")
print(f"Producer Surplus = {P_e * Q_e:.4f} - {integral_supply:.4f} = {producer_surplus:.4f}\n")

# Step 5: Calculate Total Welfare
total_welfare = consumer_surplus + producer_surplus

print("--- Total Welfare Calculation ---")
print("Total Welfare = Consumer Surplus + Producer Surplus")
print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f} = {total_welfare:.4f}")

# Final Answer
final_answer = total_welfare
print(f"\n<<<The final answer is {final_answer:.4f}>>>")