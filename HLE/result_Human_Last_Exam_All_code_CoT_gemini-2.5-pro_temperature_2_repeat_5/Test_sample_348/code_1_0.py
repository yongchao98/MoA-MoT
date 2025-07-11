import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Define the supply and demand functions
def supply(q):
    # The function is defined for q**3 - 2 > 0
    if isinstance(q, np.ndarray):
        # Handle arrays for fsolve
        # Return a large number for invalid inputs to guide the solver
        return np.where(q**3 > 2, np.log(q**3 - 2), 1e9)
    # Handle single values for quad
    if q**3 <= 2:
        return float('inf') # Integral from Q_min, so this won't be hit
    return np.log(q**3 - 2)

def demand(q):
    return 18 * np.exp(-np.arctan(q))

# Define the equation to solve for equilibrium (Demand - Supply = 0)
def equilibrium_equation(q):
    return demand(q) - supply(q)

# 1. Find the equilibrium quantity (Q_e)
# An initial guess is needed. Let's try 4.
initial_guess_q = 4.0
equilibrium_quantity, = fsolve(equilibrium_equation, initial_guess_q)

# Calculate the equilibrium price (P_e) using the demand function
equilibrium_price = demand(equilibrium_quantity)

# 2. Calculate the Total Welfare
# Minimum quantity for which supply is defined
q_min_supply = (2)**(1/3)

# Integral of the demand curve from 0 to Q_e
integral_demand, _ = quad(demand, 0, equilibrium_quantity)

# Integral of the supply curve from Q_min to Q_e
integral_supply, _ = quad(supply, q_min_supply, equilibrium_quantity)

# Total Welfare = Integral(Demand) - Integral(Supply)
total_welfare = integral_demand - integral_supply

# --- Output the results ---
print(f"Step 1: Find Market Equilibrium")
print(f"The equilibrium quantity (Qe) is: {equilibrium_quantity:.4f}")
print(f"The equilibrium price (Pe) is: {equilibrium_price:.4f}\n")

print(f"Step 2: Calculate Total Welfare")
print(f"The formula for Total Welfare is:")
print(f"TW = Integral from 0 to {equilibrium_quantity:.4f} of (18 * e^(-arctan(Q))) dQ - Integral from {q_min_supply:.4f} to {equilibrium_quantity:.4f} of (ln(Q^3 - 2)) dQ\n")

print("Solving the integrals:")
print(f"Integral of Demand from 0 to {equilibrium_quantity:.4f} = {integral_demand:.4f}")
print(f"Integral of Supply from {q_min_supply:.4f} to {equilibrium_quantity:.4f} = {integral_supply:.4f}\n")

print("Final Calculation:")
print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
print(f"Total Welfare = {total_welfare:.4f}")

<<<38.3094>>>