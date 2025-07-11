import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# This script calculates the total welfare based on given supply and demand functions.
#
# Given functions:
# Supply: P = ln(Q_S^3 - 2)
# Demand: P = 18 * e^(-arctan(Q_D))

# --- Step 1: Define the supply and demand functions ---

def demand_P(Q):
    """Calculates the price from the demand function for a given quantity Q."""
    return 18 * np.exp(-np.arctan(Q))

def supply_P(Q):
    """
    Calculates the price from the supply function for a given quantity Q.
    This function is only valid for Q > 2^(1/3).
    """
    return np.log(Q**3 - 2)

# --- Step 2: Find the market equilibrium (Q_e, P_e) ---

def equilibrium_equation(Q):
    """
    Defines the equation to solve for equilibrium.
    Equilibrium occurs when supply_P(Q) = demand_P(Q).
    The root of this function is the equilibrium quantity.
    """
    return supply_P(Q) - demand_P(Q)

# We need an initial guess to find the root. A quick analysis shows the
# intersection occurs between Q=4 and Q=5. We'll use 4.5 as a guess.
try:
    Q_e_solution = fsolve(equilibrium_equation, 4.5)
    Q_e = Q_e_solution[0]
    # Calculate the equilibrium price using the demand function
    P_e = demand_P(Q_e)

    # --- Step 3: Calculate the Total Welfare ---

    # Total Welfare = [Area under Demand curve] - [Area under Supply curve]
    # The supply curve is economically meaningful (Price >= 0) for Q >= 3^(1/3).
    # This is the starting quantity for the supply integral.
    Q_start_supply = (3.0)**(1.0/3.0)

    # Calculate the area under the demand curve from Q=0 to Q_e
    integral_demand, _ = quad(demand_P, 0, Q_e)

    # Calculate the area under the supply curve from its start (Q_start_supply) to Q_e
    integral_supply, _ = quad(supply_P, Q_start_supply, Q_e)

    # Total welfare is the difference between these two areas
    total_welfare = integral_demand - integral_supply


    # --- Step 4: Output the results ---

    print("--- Market Equilibrium ---")
    print(f"By solving ln(Q^3 - 2) = 18 * e^(-arctan(Q)), we find:")
    print(f"Equilibrium Quantity (Q_e): {Q_e:.4f}")
    print(f"Equilibrium Price (P_e):   {P_e:.4f}")
    print("\n--- Total Welfare Calculation ---")
    print("Total Welfare (TW) is calculated as the area between the demand and supply curves.")
    print("TW = (Integral of Demand from 0 to Q_e) - (Integral of Supply from its start to Q_e)")
    print(f"The starting quantity for supply (where Price=0) is 3^(1/3) ≈ {Q_start_supply:.4f}")
    print(f"\nArea under Demand = ∫ from 0 to {Q_e:.4f} [18 * e^(-arctan(Q))] dQ = {integral_demand:.4f}")
    print(f"Area under Supply = ∫ from {Q_start_supply:.4f} to {Q_e:.4f} [ln(Q^3 - 2)] dQ = {integral_supply:.4f}")
    print("\n--- Final Equation and Result ---")
    print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")

except Exception as e:
    print(f"An error occurred during calculation: {e}")
    print("Could not find the equilibrium point. Please check the functions and initial guess.")

<<<28.7188>>>