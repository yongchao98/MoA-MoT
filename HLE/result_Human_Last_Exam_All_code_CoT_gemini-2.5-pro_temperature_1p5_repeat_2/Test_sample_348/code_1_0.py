import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Define the supply and demand functions
def demand(Q):
  """Demand function: P = 18 * e^(-arctan(Q))"""
  # Add a small epsilon to Q to avoid issues with arctan(0) if necessary, although it's fine.
  return 18 * np.exp(-np.arctan(Q))

def supply(Q):
  """Supply function: P = ln(Q^3 - 2)"""
  # The function is defined for Q^3 - 2 > 0.
  # Add a check to prevent math errors for invalid Q values.
  if Q**3 - 2 <= 0:
    return -np.inf # Return a value that indicates it's outside the domain
  return np.log(Q**3 - 2)

# Define the function to find the root (equilibrium quantity)
# f(Q) = Demand(Q) - Supply(Q). We want to find Q where f(Q) = 0.
def equilibrium_equation(Q):
  return demand(Q) - supply(Q)

# --- Step 1: Find Market Equilibrium ---
# An initial guess for Q_e. Based on plotting or simple checks, Q is likely between 4 and 5.
initial_guess_Q = 4.5
# Use fsolve to find the equilibrium quantity Q_e
Q_e = fsolve(equilibrium_equation, initial_guess_Q)[0]
# Calculate the equilibrium price P_e using the demand function
P_e = demand(Q_e)

# --- Step 2: Calculate the integral of the demand function ---
# Integrate Demand(Q) from 0 to Q_e
integral_demand, _ = quad(demand, 0, Q_e)

# --- Step 3: Calculate the integral of the supply function ---
# The supply curve is defined for Q > cuberoot(2)
supply_start_Q = np.cbrt(2)
# Integrate Supply(Q) from its starting point to Q_e
integral_supply, _ = quad(supply, supply_start_Q, Q_e)

# --- Step 4: Calculate Total Welfare ---
# Total Welfare = (Integral of Demand) - (Integral of Supply)
total_welfare = integral_demand - integral_supply

# --- Output the results ---
print(f"Market Equilibrium:")
print(f"Equilibrium Quantity (Q_e): {Q_e:.4f}")
print(f"Equilibrium Price (P_e): {P_e:.4f}")
print("\nCalculating Total Welfare:")
print(f"Total Welfare = (Integral of Demand from 0 to {Q_e:.4f}) - (Integral of Supply from {supply_start_Q:.4f} to {Q_e:.4f})")
print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
print(f"Total Welfare = {total_welfare:.4f}")

# Final Answer
print(f"\nFinal Answer: The total welfare is {total_welfare:.4f}")
<<<45.5088>>>