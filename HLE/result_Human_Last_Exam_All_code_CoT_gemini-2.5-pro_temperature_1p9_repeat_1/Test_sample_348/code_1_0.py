import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad

# 1. Define supply and demand functions
def supply_func(q):
    """Defines the supply function P = ln(Q^3 - 2)"""
    if q**3 <= 2:
        return np.nan # Not defined for q <= 2^(1/3)
    return np.log(q**3 - 2)

def demand_func(q):
    """Defines the demand function P = 18 * e^(-arctan(Q))"""
    return 18 * np.exp(-np.arctan(q))

# 2. Find market equilibrium
# We need to solve demand_func(q) - supply_func(q) = 0
def equilibrium_difference(q):
    # This check ensures the function is well-behaved for the solver
    if q**3 <= 2:
        # Return a large positive value if q is too low, as demand price will be higher than supply price
        return 1.0 
    return demand_func(q) - supply_func(q)

# Solving numerically. A quick check shows the root is between 4 and 5.
# equilibrium_difference(4) is positive, equilibrium_difference(5) is negative.
solution = root_scalar(equilibrium_difference, bracket=[4, 5], method='brentq')
Q_e = solution.root
P_e = demand_func(Q_e)

# 3. Calculate Consumer Surplus (CS)
# CS = Integral(Demand) from 0 to Q_e - P_e * Q_e
integral_demand, _ = quad(demand_func, 0, Q_e)
consumer_surplus = integral_demand - P_e * Q_e

# 4. Calculate Producer Surplus (PS)
# PS = P_e * Q_e - Integral(Supply) from Q_min to Q_e
# The supply function is defined for Q > 2^(1/3)
q_min_supply = (2)**(1/3)
integral_supply, _ = quad(supply_func, q_min_supply, Q_e)
producer_surplus = P_e * Q_e - integral_supply

# 5. Calculate Total Welfare (TW)
total_welfare = consumer_surplus + producer_surplus

# Print all the calculated values
print(f"Market Equilibrium:")
print(f"Equilibrium Quantity (Qe): {Q_e:.4f}")
print(f"Equilibrium Price (Pe): {P_e:.4f}\n")

print(f"Surplus Calculation:")
print(f"Consumer Surplus: {consumer_surplus:.4f}")
print(f"Producer Surplus: {producer_surplus:.4f}\n")

print("The Final Equation for Total Welfare is:")
print(f"Total Welfare = Consumer Surplus + Producer Surplus")
print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f}")
print(f"Total Welfare = {total_welfare:.4f}")

<<<22.0122>>>