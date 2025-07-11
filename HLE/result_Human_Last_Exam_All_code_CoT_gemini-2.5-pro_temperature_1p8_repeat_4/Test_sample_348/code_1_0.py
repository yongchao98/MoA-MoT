import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare for a given set of supply and demand functions.
    """
    
    # 1. Define the supply and demand functions
    def demand(q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(q))

    def supply(q):
        """Supply function: P = ln(Q^3 - 2)"""
        # Ensure the argument of log is positive to avoid math domain errors.
        if q**3 - 2 <= 0:
            return -np.inf  # Return a value that signals an invalid domain
        return np.log(q**3 - 2)

    # 2. Find the market equilibrium point (Qe, Pe)
    def equilibrium_equation(q):
        """Equation to solve for equilibrium: Demand(Q) - Supply(Q) = 0"""
        return demand(q) - supply(q)

    # The supply function is defined for Q > 2^(1/3) ~= 1.26
    # We look for a root in a bracket where the function values have opposite signs.
    # At Q=4, demand > supply. At Q=5, supply > demand. The root is between 4 and 5.
    try:
        sol = root_scalar(equilibrium_equation, bracket=[4, 5])
        if not sol.converged:
            print("Error: The equilibrium quantity could not be found.")
            return
        q_equilibrium = sol.root
        p_equilibrium = demand(q_equilibrium)
    except ValueError:
        print("Error: The provided bracket for root finding is invalid.")
        return
        
    # 3. Calculate Total Welfare
    # Total Welfare = Integral(Demand) - Integral(Supply) over their respective domains.
    
    # Integral of the demand curve from Q=0 to Qe
    integral_demand, _ = quad(demand, 0, q_equilibrium)
    
    # The lower bound for the supply curve is where Q^3 - 2 > 0 => Q > 2^(1/3)
    q_min_supply = 2**(1/3)
    
    # Integral of the supply curve from Q_min to Qe
    integral_supply, _ = quad(supply, q_min_supply, q_equilibrium)

    # Calculate total welfare
    total_welfare = integral_demand - integral_supply

    # 4. Print the results
    print(f"Market Equilibrium:")
    print(f"  - Equilibrium Quantity (Qe): {q_equilibrium:.4f}")
    print(f"  - Equilibrium Price (Pe): {p_equilibrium:.4f}\n")
    
    print("Total Welfare Calculation:")
    print(f"  - Area under Demand Curve (Consumer Benefit): {integral_demand:.4f}")
    print(f"  - Area under Supply Curve (Producer Cost): {integral_supply:.4f}\n")

    print("The final equation for Total Welfare is:")
    print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f} = {total_welfare:.4f}")

    # Return the final numerical answer in the required format
    # The string conversion is just to add the required brackets
    final_answer = f"<<<{total_welfare:.4f}>>>"
    #print(final_answer) # This would print the <<<answer>>> line

if __name__ == "__main__":
    main()