import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Define the supply and demand functions
# Supply: P = ln(Q^3 - 2)
def supply_func(q):
    # This function is only defined for Q > 2^(1/3)
    # Return a very large number for invalid inputs to help the solver
    if np.any(q <= 2**(1/3)):
        return 1e9
    return np.log(q**3 - 2)

# Demand: P = 18 * e^(-arctan(Q))
def demand_func(q):
    return 18 * np.exp(-np.arctan(q))

# Define the equation to find equilibrium by finding the root of f(Q) = 0
# f(Q) = Supply(Q) - Demand(Q)
def equilibrium_equation(q):
    return supply_func(q) - demand_func(q)

# --- Step 1: Find Market Equilibrium ---
# Use a numerical solver to find the equilibrium quantity (Q_E)
# An initial guess of Q=2.0 is reasonable
q_initial_guess = 2.0
# fsolve returns an array, so we take the first element
q_equilibrium = fsolve(equilibrium_equation, q_initial_guess)[0]

# Calculate the equilibrium price (P_E) using the demand function
p_equilibrium = demand_func(q_equilibrium)

print(f"Market Equilibrium:")
print(f"  Equilibrium Quantity (Q_E): {q_equilibrium:.4f}")
print(f"  Equilibrium Price (P_E): {p_equilibrium:.4f}\n")


# --- Step 2 & 3: Calculate Consumer and Producer Surplus ---
# The supply curve starts where Q^3 - 2 > 0, so Q > 2^(1/3)
supply_integration_start = 2**(1/3)

# Calculate the integral of the demand curve from 0 to Q_E
# quad returns a tuple (result, error), so we take the first element
integral_demand, _ = quad(demand_func, 0, q_equilibrium)

# Calculate the integral of the supply curve from its start to Q_E
integral_supply, _ = quad(supply_func, supply_integration_start, q_equilibrium)

# Calculate Consumer Surplus (CS)
consumer_surplus = integral_demand - (p_equilibrium * q_equilibrium)

# Calculate Producer Surplus (PS)
producer_surplus = (p_equilibrium * q_equilibrium) - integral_supply

print("Surplus Calculation:")
print(f"  Consumer Surplus (CS) = {consumer_surplus:.4f}")
print(f"  Producer Surplus (PS) = {producer_surplus:.4f}\n")


# --- Step 4: Calculate Total Welfare ---
# Total Welfare is the sum of consumer and producer surplus
total_welfare = consumer_surplus + producer_surplus

print("Final Calculation:")
print(f"Total Welfare = Consumer Surplus + Producer Surplus")
print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f}")
print(f"Total Welfare = {total_welfare:.4f}")

<<<22.3023>>>