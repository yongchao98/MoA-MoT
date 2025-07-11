import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad

# --- 1. Define Supply and Demand Functions ---

def demand_func(q):
    """Demand function: P = 18 * e^(-arctan(Q))"""
    return 18 * np.exp(-np.arctan(q))

def supply_func(q):
    """Supply function: P = ln(Q^3 - 2)"""
    # The function is defined only for Q > cbrt(2)
    if q**3 - 2 <= 0:
        return -np.inf  # Return a value that indicates it's outside the domain
    return np.log(q**3 - 2)

# --- 2. Find Market Equilibrium (Qe, Pe) ---

# To find equilibrium, we find the root of the difference between demand and supply
def equilibrium_equation(q):
    return demand_func(q) - supply_func(q)

# The supply function domain requires Q > ∛2 ≈ 1.26.
# By plotting or testing values, we can find a bracket for the root, e.g., [4, 5].
try:
    solution = root_scalar(equilibrium_equation, bracket=[4, 5])
    Q_e = solution.root
    # Calculate Pe using the demand function (or supply, they are equal at Qe)
    P_e = demand_func(Q_e)
    print(f"Market Equilibrium:")
    print(f"Equilibrium Quantity (Q_e): {Q_e:.4f}")
    print(f"Equilibrium Price (P_e): {P_e:.4f}\n")

    # --- 3. Calculate Consumer Surplus (CS) ---
    # CS = ∫[0 to Qe] D(q)dq - Pe * Qe
    integral_demand, _ = quad(demand_func, 0, Q_e)
    consumer_surplus = integral_demand - P_e * Q_e
    print(f"Surplus Calculation:")
    print(f"Consumer Surplus (CS): {consumer_surplus:.4f}")

    # --- 4. Calculate Producer Surplus (PS) ---
    # PS = Pe * Qe - ∫[Q_min to Qe] S(q)dq
    # The lower bound of integration for supply is the quantity where P=0, but here it's where the function is defined.
    q_min_supply = np.cbrt(2)
    integral_supply, _ = quad(supply_func, q_min_supply, Q_e)
    producer_surplus = P_e * Q_e - integral_supply
    print(f"Producer Surplus (PS): {producer_surplus:.4f}\n")

    # --- 5. Calculate and Display Total Welfare ---
    total_welfare = consumer_surplus + producer_surplus
    print(f"Total Welfare Equation (CS + PS = TW):")
    # Final output showing the numbers in the equation
    print(f"{consumer_surplus:.4f} + {producer_surplus:.4f} = {total_welfare:.4f}")

except ValueError as e:
    print(f"An error occurred during calculation: {e}")
    print("Could not find the equilibrium point in the specified bracket.")
