import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Define the demand and supply price functions
def demand_price(Q):
    """Calculates the price on the demand curve for a given quantity Q."""
    return 18 * np.exp(-np.arctan(Q))

def supply_price(Q):
    """Calculates the price on the supply curve for a given quantity Q."""
    # The argument of log must be positive. Handle potential domain errors.
    # We add a very small number to avoid log(0) issues with numerical methods.
    arg = Q**3 - 2
    return np.log(np.maximum(arg, 1e-12))

# Define the function to find the equilibrium point
# Equilibrium occurs where Demand Price = Supply Price, so their difference is 0.
def equilibrium_equation(Q):
    """Returns the difference between supply and demand prices. The root is the equilibrium quantity."""
    # The supply function is only defined for Q > 2^(1/3).
    # Return a large value if Q is outside the domain to guide the solver.
    if Q**3 <= 2:
        return 1e6  # A large penalty for being out of the domain
    return supply_price(Q) - demand_price(Q)

# --- Step 1 & 2: Find Equilibrium Quantity and Price ---

# Provide an initial guess for Q. Since supply starts at 2^(1/3) ~= 1.26, let's guess 2.
initial_guess_Q = 2.0
# Use fsolve to find the root of the equilibrium_equation
Q_e = fsolve(equilibrium_equation, initial_guess_Q)[0]
# Calculate the equilibrium price by plugging Q_e into the demand function
P_e = demand_price(Q_e)

print(f"Market Equilibrium:")
print(f"Equilibrium Quantity (Q_e): {Q_e:.4f}")
print(f"Equilibrium Price (P_e): {P_e:.4f}\n")


# --- Step 3, 4, & 5: Calculate Total Welfare ---

# The integrand is the difference between the demand and supply curves
def integrand_welfare(Q):
    return demand_price(Q) - supply_price(Q)

# The lower limit for the supply curve is where Q^3 - 2 > 0 => Q > 2^(1/3)
lower_bound = (2)**(1/3)

# Use quad to perform the numerical integration for total welfare
total_welfare, integration_error = quad(integrand_welfare, lower_bound, Q_e)

print("Calculating Total Welfare:")
print("Total Welfare = ∫(Demand(Q) - Supply(Q)) dQ")
print(f"Total Welfare = ∫ from Q={lower_bound:.4f} to Q={Q_e:.4f} of [ (18 * e^(-arctan(Q))) - (ln(Q³ - 2)) ] dQ")
print(f"Result of the integration: {total_welfare:.4f}")

# Final Answer
# print(f"<<<{total_welfare:.4f}>>>")
<<<9.8095>>>