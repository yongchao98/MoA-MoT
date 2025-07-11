import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Define the supply and demand functions
def supply_func(Q):
    # The function is defined for Q^3 - 2 > 0
    if np.any(Q**3 - 2 <= 0):
        return np.inf  # Return a large number for invalid inputs
    return np.log(Q**3 - 2)

def demand_func(Q):
    return 18 * np.exp(-np.arctan(Q))

# Define the function to find the root of, which represents equilibrium
def equilibrium_equation(Q):
    return demand_func(Q) - supply_func(Q)

# 1. Find market equilibrium
# Initial guess for the solver, based on preliminary analysis
initial_guess_Q = 4.5
# Use fsolve to find the equilibrium quantity Q*
Q_equilibrium = fsolve(equilibrium_equation, initial_guess_Q)[0]

# Calculate the equilibrium price P*
P_equilibrium = demand_func(Q_equilibrium)

# 2. Calculate Total Welfare
# The supply curve starts producing at a positive price when ln(Q^3-2) >= 0 => Q^3-2 >= 1 => Q^3 >= 3
Q_supply_start = 3**(1/3)

# Define the integrands for scipy.integrate.quad
integrand_demand = demand_func
integrand_supply = supply_func

# Calculate the integral of the demand curve from 0 to Q_equilibrium
integral_demand, _ = quad(integrand_demand, 0, Q_equilibrium)

# Calculate the integral of the supply curve from where it starts to Q_equilibrium
integral_supply, _ = quad(integrand_supply, Q_supply_start, Q_equilibrium)

# Total welfare is the area between the curves
total_welfare = integral_demand - integral_supply

# 3. Print the results
print(f"Market Equilibrium:")
print(f"Equilibrium Quantity (Q*): {Q_equilibrium}")
print(f"Equilibrium Price (P*): {P_equilibrium}\n")

print("Total Welfare Calculation:")
# Print the final equation with the computed values
print("Total Welfare = ∫[0 to Q*] (18 * e^(-arctan(Q))) dQ - ∫[Q_supply_start to Q*] (ln(Q^3 - 2)) dQ")
print(f"Total Welfare = ∫[0 to {Q_equilibrium:.4f}] (18 * e^(-arctan(Q))) dQ - ∫[{Q_supply_start:.4f} to {Q_equilibrium:.4f}] (ln(Q^3 - 2)) dQ\n")

print(f"Final calculated Total Welfare: {total_welfare}")

print(f"\n<<<Total Welfare = {total_welfare:.4f}>>>")
<<<Total Welfare = 16.9691>>>