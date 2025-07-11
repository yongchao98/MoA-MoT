import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Step 1: Define the supply and demand functions
# Supply: P = ln(Q^3 - 2)
def supply_price(Q):
    # The function is defined for Q^3 - 2 > 0.
    # We need to handle cases where the argument to log is non-positive for the solver.
    if np.any(Q**3 - 2 <= 0):
        # Return a large value to push the solver away from this invalid region
        return 1e9
    return np.log(Q**3 - 2)

# Demand: P = 18 * e^(-arctan(Q))
def demand_price(Q):
    return 18 * np.exp(-np.arctan(Q))

# Step 2: Find the market equilibrium
# Equilibrium occurs when supply price equals demand price.
# We need to find the root of the equation: demand_price(Q) - supply_price(Q) = 0
def equilibrium_equation(Q):
    return demand_price(Q) - supply_price(Q)

# Provide an initial guess for the solver. A quick check shows the
# equilibrium quantity is likely between 3 and 4.
initial_guess = 3.5
# fsolve returns an array, so we take the first element
Q_e = fsolve(equilibrium_equation, initial_guess)[0]

# Step 3: Calculate the components for the total welfare equation
# The integral of the demand curve is from Q=0 to Q_e
integral_demand, _ = quad(demand_price, 0, Q_e)

# The supply curve is only defined for Q > 2^(1/3).
# So, the integral of the supply curve is from Q_start = 2^(1/3) to Q=Q_e.
Q_supply_start = np.cbrt(2)
integral_supply, _ = quad(supply_price, Q_supply_start, Q_e)

# Step 4: Calculate the final total welfare
total_welfare = integral_demand - integral_supply

# Step 5: Print the final equation with all its components
print("The total welfare is calculated as the area between the demand and supply curves.")
print("Total Welfare = (Integral of Demand) - (Integral of Supply)")
print("\nHere are the numbers for the final equation:")
print(f"Integral of Demand from 0 to Q_e={Q_e:.4f}: {integral_demand:.4f}")
print(f"Integral of Supply from Q_start={Q_supply_start:.4f} to Q_e={Q_e:.4f}: {integral_supply:.4f}")
print("\nFinal Equation:")
print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f} = {total_welfare:.4f}")

<<<21.0399>>>