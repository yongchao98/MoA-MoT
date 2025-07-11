import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Define the supply and demand functions
def demand_p(q):
    """Calculates price based on demand quantity."""
    return 18 * np.exp(-np.arctan(q))

def supply_p(q):
    """Calculates price based on supply quantity."""
    # Ensure the input is valid for the natural logarithm
    q_safe = np.asarray(q)
    if np.any(q_safe**3 - 2 <= 0):
        # Return a large number for out-of-domain values to guide the solver
        return np.inf 
    return np.log(q**3 - 2)

# Define the equation to find the equilibrium point (Supply Price - Demand Price = 0)
def equilibrium_equation(q):
    """Difference between supply and demand prices. The root is the equilibrium quantity."""
    return supply_p(q) - demand_p(q)

# --- Step 1: Find Market Equilibrium ---
# The supply function is defined for Q^3 > 2, so Q > 2^(1/3) ~= 1.26
# We need an initial guess for the solver in the valid domain.
initial_guess = 2.0 
q_equilibrium_array = fsolve(equilibrium_equation, initial_guess)
q_equilibrium = q_equilibrium_array[0]
p_equilibrium = demand_p(q_equilibrium)

print("--- Step 1: Find Market Equilibrium ---")
print("We solve the equation: P_Supply = P_Demand")
print("ln(Q^3 - 2) = 18 * e^(-arctan(Q))")
print(f"Using a numerical solver, we find the equilibrium point:")
print(f"Equilibrium Quantity (Q_E): {q_equilibrium:.4f}")
print(f"Equilibrium Price (P_E): {p_equilibrium:.4f}\n")

# --- Step 2: Calculate Total Welfare ---
# Minimum quantity where the supply curve is defined
q_min_supply = (2)**(1/3)

# Integrate the demand function from 0 to Q_E
integral_demand, _ = quad(demand_p, 0, q_equilibrium)

# Integrate the supply function from its starting point (Q_min_supply) to Q_E
integral_supply, _ = quad(supply_p, q_min_supply, q_equilibrium)

# Total welfare is the area under the demand curve minus the area under the supply curve
total_welfare = integral_demand - integral_supply

print("--- Step 2: Calculate Total Welfare ---")
print("Total Welfare (TW) is calculated as the area between the demand and supply curves.")
print("TW = (Area under Demand Curve from 0 to Q_E) - (Area under Supply Curve from Q_min to Q_E)")
print("The final equation is:")
print(f"TW = ∫[from 0 to {q_equilibrium:.4f}] (18 * e^(-arctan(Q))) dQ - ∫[from {q_min_supply:.4f} to {q_equilibrium:.4f}] (ln(Q^3 - 2)) dQ")
print("\nCalculating the values of the integrals:")
print(f"Value of Demand Integral: {integral_demand:.4f}")
print(f"Value of Supply Integral: {integral_supply:.4f}")
print("\n--- Final Calculation ---")
print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
print(f"Total Welfare = {total_welfare:.4f}")

<<<38.7904>>>