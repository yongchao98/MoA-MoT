import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def solve_and_print_welfare():
    """
    This function calculates and prints the steps to find the total welfare
    given supply and demand functions.
    """
    # Step 1: Define the supply and demand functions
    def supply_p(q):
        # The function is defined for q^3 - 2 > 0.
        # Handle invalid inputs for the solver and integrator.
        if isinstance(q, (np.ndarray, list)):
            q_safe = np.asarray(q)
            # Create a boolean mask for valid domain
            valid_mask = q_safe**3 > 2
            result = np.full_like(q_safe, np.nan)
            # Calculate only for values within the domain
            result[valid_mask] = np.log(q_safe[valid_mask]**3 - 2)
            return result
        else: # Handle single float/int
            if q**3 > 2:
                return np.log(q**3 - 2)
            else:
                return np.nan

    def demand_p(q):
        return 18 * np.exp(-np.arctan(q))

    # Step 2: Find the market equilibrium point (Q_e, P_e)
    # We solve Supply(Q) = Demand(Q), or Supply(Q) - Demand(Q) = 0
    def equilibrium_equation(q):
        s_price = supply_p(q)
        # If solver tries a value outside the domain, s_price will be NaN.
        # Return a large number to guide the solver away from this region.
        if np.isnan(s_price):
            return 1e10
        return s_price - demand_p(q)

    # Based on plotting or testing values, a good initial guess is between 4 and 5.
    initial_guess = 4.5
    q_equilibrium = fsolve(equilibrium_equation, initial_guess)[0]
    p_equilibrium = demand_p(q_equilibrium)

    # Step 3: Calculate surpluses
    # The supply function's domain starts at Q_min where Q^3 - 2 = 0
    q_min_supply = (2)**(1/3)

    # Consumer Surplus (CS) = integral from 0 to Q_e of (Demand(Q) - P_e) dQ
    cs_integrand = lambda q: demand_p(q) - p_equilibrium
    consumer_surplus, _ = quad(cs_integrand, 0, q_equilibrium)

    # Producer Surplus (PS) = integral from Q_min_supply to Q_e of (P_e - Supply(Q)) dQ
    ps_integrand = lambda q: p_equilibrium - supply_p(q)
    producer_surplus, _ = quad(ps_integrand, q_min_supply, q_equilibrium)

    # Total Welfare (TW) = CS + PS
    total_welfare = consumer_surplus + producer_surplus

    # Step 4: Print the results and the "equations" with numbers
    print("To find the total welfare, we follow these steps:")
    print("\n1. Find the market equilibrium where Supply equals Demand.")
    print(f"   Supply: P = ln(Q³ - 2)")
    print(f"   Demand: P = 18 * e^(-arctan(Q))")
    print(f"   Solving ln(Q³ - 2) = 18 * e^(-arctan(Q)) numerically gives:")
    print(f"   Equilibrium Quantity (Q_e) = {q_equilibrium:.4f}")
    print(f"   Equilibrium Price (P_e) = {p_equilibrium:.4f}")

    print("\n2. Calculate the Consumer Surplus (CS).")
    print("   CS = ∫[from 0 to Q_e] (Demand(Q) - P_e) dQ")
    print(f"   CS = ∫[from 0 to {q_equilibrium:.4f}] (18 * e^(-arctan(Q)) - {p_equilibrium:.4f}) dQ")
    print(f"   Calculated Consumer Surplus (CS) ≈ {consumer_surplus:.4f}")

    print("\n3. Calculate the Producer Surplus (PS).")
    print(f"   The supply curve starts at Q_min = 2^(1/3) ≈ {q_min_supply:.4f}")
    print("   PS = ∫[from Q_min to Q_e] (P_e - Supply(Q)) dQ")
    print(f"   PS = ∫[from {q_min_supply:.4f} to {q_equilibrium:.4f}] ({p_equilibrium:.4f} - ln(Q³ - 2)) dQ")
    print(f"   Calculated Producer Surplus (PS) ≈ {producer_surplus:.4f}")

    print("\n4. Calculate the Total Welfare (TW).")
    print("   TW = CS + PS")
    print(f"   TW = {consumer_surplus:.4f} + {producer_surplus:.4f}")
    print(f"   Final Total Welfare (TW) ≈ {total_welfare:.4f}")

    # The final answer to be extracted
    # Do not remove the next line
    print(f"\n<<<{total_welfare:.4f}>>>")

if __name__ == '__main__':
    solve_and_print_welfare()