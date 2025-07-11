import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Step 1: Define the supply and demand functions
def demand(q):
    """Demand function: P = 18 * e^(-arctan(Q))"""
    return 18 * np.exp(-np.arctan(q))

def supply(q):
    """Supply function: P = ln(Q^3 - 2)"""
    # The function is defined for Q^3 - 2 > 0, which means Q > cuberoot(2)
    # To handle integration, we return a large number for values outside the domain
    # This prevents errors but fsolve should find a root within the valid domain.
    if isinstance(q, (np.ndarray, list)):
        q_safe = np.array(q)
        q_safe[q_safe**3 <= 2] = np.inf
        return np.log(q_safe**3 - 2)
    else:
        if q**3 <= 2:
            return np.inf # Not defined
        return np.log(q**3 - 2)

# Step 2: Find the market equilibrium quantity (Q_E)
# We need to solve the equation: demand(q) - supply(q) = 0
def equilibrium_equation(q):
    """Difference between demand and supply functions"""
    return demand(q) - supply(q)

# An initial guess is needed for the solver. Let's try Q=4.
# demand(4) = 4.78, supply(4) = 4.12. They are close.
initial_guess_q = 4.0
# Use fsolve to find the root, which is the equilibrium quantity
equilibrium_quantity = fsolve(equilibrium_equation, initial_guess_q)[0]

# Step 3: Calculate the equilibrium price (P_E)
# Substitute Q_E back into the demand function
equilibrium_price = demand(equilibrium_quantity)

# Step 4: Calculate Total Welfare
# Total Welfare = Integral(Demand) from 0 to Q_E - Integral(Supply) from its start to Q_E
# The supply curve starts where Q^3 - 2 > 0 => Q > cuberoot(2)
supply_start_q = np.cbrt(2)

# Calculate the integral of the demand curve from 0 to Q_E
integral_demand, _ = quad(demand, 0, equilibrium_quantity)

# Calculate the integral of the supply curve from its starting quantity to Q_E
integral_supply, _ = quad(supply, supply_start_q, equilibrium_quantity)

# Calculate the final total welfare
total_welfare = integral_demand - integral_supply

# Step 5: Output the results
print(f"Market Equilibrium:")
print(f"  Equilibrium Quantity (Q_E) = {equilibrium_quantity:.4f}")
print(f"  Equilibrium Price (P_E) = {equilibrium_price:.4f}\n")
print(f"Total Welfare Calculation:")
print(f"The final equation for Total Welfare is: Integral(Demand) - Integral(Supply)")
print(f"  Integral of Demand from 0 to {equilibrium_quantity:.4f} = {integral_demand:.4f}")
print(f"  Integral of Supply from {supply_start_q:.4f} to {equilibrium_quantity:.4f} = {integral_supply:.4f}")
print(f"  Total Welfare = {integral_demand:.4f} - {integral_supply:.4f} = {total_welfare:.4f}\n")
print(f"The total welfare is: {total_welfare:.4f}")

<<<11.5303>>>