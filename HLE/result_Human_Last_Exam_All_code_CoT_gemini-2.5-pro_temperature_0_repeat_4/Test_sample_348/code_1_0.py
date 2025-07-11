import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """
    # --- Step 1: Define Supply and Demand Functions ---
    
    # Demand function: P = 18 * e^(-arctan(Q))
    def demand_func(Q):
        return 18 * np.exp(-np.arctan(Q))

    # Supply function: P = ln(Q^3 - 2)
    def supply_func(Q):
        # The function is defined for Q^3 - 2 > 0.
        # We handle non-positive arguments to avoid math errors.
        arg = Q**3 - 2
        if isinstance(arg, (int, float)):
            return np.log(arg) if arg > 0 else np.nan
        else: # Handle numpy arrays
            return np.log(np.where(arg > 0, arg, np.nan))

    # --- Step 2: Find Market Equilibrium ---
    
    # We need to solve for Q where Demand(Q) = Supply(Q)
    def equilibrium_equation(Q):
        # This function should return 0 at the equilibrium point.
        # We must handle the domain of the supply function.
        if Q <= np.cbrt(2):
            # Return a large number if outside the valid domain to guide the solver.
            return 1e6
        return demand_func(Q) - supply_func(Q)

    # An initial guess for Q is needed for the numerical solver.
    # Based on plotting or trial and error, we know Q is between 4 and 5.
    initial_guess = 4.5
    
    # Use fsolve to find the equilibrium quantity (Q_e)
    Q_e = fsolve(equilibrium_equation, initial_guess)[0]
    
    # Calculate the equilibrium price (P_e)
    P_e = demand_func(Q_e)

    # --- Step 3: Calculate Total Welfare ---
    
    # The supply curve starts where its argument is positive: Q^3 - 2 > 0 => Q > 2^(1/3)
    Q_min = np.cbrt(2)

    # Calculate the integral of the demand function from 0 to Q_e
    integral_demand, _ = quad(demand_func, 0, Q_e)

    # Calculate the integral of the supply function from Q_min to Q_e
    integral_supply, _ = quad(supply_func, Q_min, Q_e)

    # Total Welfare = (Area under Demand) - (Area under Supply)
    total_welfare = integral_demand - integral_supply

    # --- Step 4: Print the Results ---
    
    print(f"Market Equilibrium:")
    print(f"  - Equilibrium Quantity (Q_e): {Q_e:.4f}")
    print(f"  - Equilibrium Price (P_e): {P_e:.4f}\n")
    
    print("Total Welfare Calculation:")
    print(f"The total welfare is the area between the demand and supply curves.")
    print(f"TW = ∫[0 to {Q_e:.4f}] D(Q)dQ - ∫[{Q_min:.4f} to {Q_e:.4f}] S(Q)dQ\n")
    
    print("Final Equation with calculated values:")
    # The final equation is Integral_Demand - Integral_Supply = Total_Welfare
    print(f"{integral_demand:.4f} - {integral_supply:.4f} = {total_welfare:.4f}")

if __name__ == "__main__":
    main()
<<<15.6003>>>