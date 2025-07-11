import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Step 1: Define the supply and demand functions
# Demand: P = 18 * exp(-arctan(Q))
def demand_func(q):
    return 18 * np.exp(-np.arctan(q))

# Supply: P = ln(Q^3 - 2)
def supply_func(q):
    # This function is defined only for Q > 2^(1/3)
    # np.log handles arrays, but we ensure the input is valid for the problem context
    # In integration, quad handles the domain, but fsolve needs guidance
    with np.errstate(invalid='ignore'): # Suppress log(negative) warning during search
        result = np.log(q**3 - 2)
    return result

# Create a function representing Demand - Supply = 0 to find the equilibrium point
def equilibrium_equation(q):
    return demand_func(q) - supply_func(q)

# Find the equilibrium quantity (Q_e) using a numerical solver
# An initial guess of 4.5 is chosen based on preliminary analysis
initial_guess_q = 4.5
Q_e = fsolve(equilibrium_equation, initial_guess_q)[0]

# Calculate the equilibrium price (P_e) using the demand function
P_e = demand_func(Q_e)

# Step 2: Calculate Consumer Surplus (CS)
# CS = Integral from 0 to Q_e of (Demand(Q)) dQ - (P_e * Q_e)
consumer_integral, _ = quad(demand_func, 0, Q_e)
consumer_surplus = consumer_integral - (P_e * Q_e)

# Step 3: Calculate Producer Surplus (PS)
# Find the minimum quantity where the supply function is defined (Q^3 - 2 > 0)
Q_min = np.cbrt(2)
# PS = (P_e * Q_e) - Integral from Q_min to Q_e of (Supply(Q)) dQ
producer_integral, _ = quad(supply_func, Q_min, Q_e)
producer_surplus = (P_e * Q_e) - producer_integral

# Step 4: Calculate Total Welfare
# Total Welfare = Consumer Surplus + Producer Surplus
total_welfare = consumer_surplus + producer_surplus

# Print the results
print(f"Market Equilibrium:")
print(f"  Equilibrium Quantity (Q_e): {Q_e:.4f}")
print(f"  Equilibrium Price (P_e):    {P_e:.4f}\n")

print(f"Calculating Total Welfare:")
print(f"  Consumer Surplus (CS) = {consumer_surplus:.4f}")
print(f"  Producer Surplus (PS) = {producer_surplus:.4f}")
print(f"  Total Welfare (CS + PS) = {total_welfare:.4f}")

# The line below is for the final answer detection
# <<< {total_welfare:.4f} >>>
<<<45.8569>>>