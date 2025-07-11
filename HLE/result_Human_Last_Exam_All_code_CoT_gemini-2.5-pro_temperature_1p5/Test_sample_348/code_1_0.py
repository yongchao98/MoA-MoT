import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# 1. Define the supply and demand functions
def demand_P(Q):
    """Demand function P = 18 * e^(-arctan(Q))"""
    return 18 * np.exp(-np.arctan(Q))

def supply_P(Q):
    """Supply function P = ln(Q^3 - 2)"""
    # Use np.maximum to avoid domain errors during solving
    # The actual integration will use the correct lower bound.
    return np.log(np.maximum(Q**3 - 2, 1e-9))

# 2. Find market equilibrium (Qe, Pe)
# Define the equation to solve: supply_P(Q) - demand_P(Q) = 0
def equilibrium_equation(Q):
    return supply_P(Q) - demand_P(Q)

# Initial guess for Q. From analysis, we know Q > 1.26. Let's guess around 4.
initial_guess_Q = 4.0
# Use fsolve to find the root, which is the equilibrium quantity Qe
Q_e_solution = fsolve(equilibrium_equation, initial_guess_Q)
Qe = Q_e_solution[0]

# Calculate the equilibrium price Pe
Pe = demand_P(Qe)

# 3. Determine the integration limit for the supply curve
# Supply is defined for Q^3 - 2 > 0 => Q > 2^(1/3)
Q_min = 2**(1/3)

# 4. Calculate the components of Total Welfare
# Integral of the demand curve from 0 to Qe
integral_demand, _ = quad(demand_P, 0, Qe)

# Integral of the supply curve from Q_min to Qe
integral_supply, _ = quad(supply_P, Q_min, Qe)

# Calculate total welfare
total_welfare = integral_demand - integral_supply

# 5. Output the results
print("--- Market Equilibrium ---")
print(f"Equilibrium Quantity (Qe): {Qe:.4f}")
print(f"Equilibrium Price (Pe): {Pe:.4f}")
print("\n--- Total Welfare Calculation ---")
print(f"The lower bound for the supply curve (Q_min = 2^(1/3)) is: {Q_min:.4f}")
print("Total Welfare (TW) is calculated as: ∫[0 to Qe] P_D(Q)dQ - ∫[Q_min to Qe] P_S(Q)dQ\n")
print("Equation with calculated values:")
print(f"Integral of Demand from 0 to {Qe:.4f} = {integral_demand:.4f}")
print(f"Integral of Supply from {Q_min:.4f} to {Qe:.4f} = {integral_supply:.4f}")
print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
print(f"Total Welfare = {total_welfare:.4f}")
print(f'<<<{total_welfare:.4f}>>>')