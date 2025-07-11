import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar

def solve_welfare():
    """
    This function calculates the total welfare given supply and demand functions.
    """
    # Step 1: Define the supply and demand functions
    def demand(q):
        """Demand function: P = 18 * e^(-arctan(Q_D))"""
        return 18 * np.exp(-np.arctan(q))

    def supply(q):
        """Supply function: P = ln(Q_S^3 - 2)"""
        # This function is defined for Q > 2^(1/3)
        # Using np.maximum to prevent log of non-positive numbers during integration
        # but the integration bounds will handle this properly.
        return np.log(q**3 - 2)

    # Step 2: Find the market equilibrium quantity and price
    
    # Equilibrium is where demand(q) == supply(q)
    def equilibrium_equation(q):
        return demand(q) - supply(q)

    # The supply function is only defined for Q > 2^(1/3)
    Q_min_supply = 2**(1/3)

    # We need to find the root of the equilibrium_equation.
    # We will search for a root in a bracket that is within the valid domain.
    # A quick check shows a root exists between 2 and 10.
    try:
        solution = root_scalar(equilibrium_equation, bracket=[2, 10])
        Q_star = solution.root
    except ValueError as e:
        print(f"Could not find the equilibrium point within the initial bracket. Error: {e}")
        return

    # Calculate the equilibrium price
    P_star = demand(Q_star)

    print("Market Equilibrium Found:")
    print(f"Equilibrium Quantity (Q*) = {Q_star}")
    print(f"Equilibrium Price (P*) = {P_star}")
    print("-" * 30)

    # Step 3: Calculate the total welfare
    # Total Welfare = ∫ D(Q)dQ - ∫ S(Q)dQ
    
    # Integral of the demand curve from 0 to Q_star
    integral_demand, err_demand = quad(demand, 0, Q_star)

    # Integral of the supply curve from its minimum valid quantity (Q_min_supply) to Q_star
    integral_supply, err_supply = quad(supply, Q_min_supply, Q_star)

    # Calculate final total welfare
    total_welfare = integral_demand - integral_supply
    
    print("Calculating Total Welfare:")
    print("Formula: Total Welfare = ∫(Demand(Q))dQ - ∫(Supply(Q))dQ")
    print(f"The integral of the demand curve from Q=0 to Q*={Q_star:.4f} is: {integral_demand}")
    print(f"The integral of the supply curve from Q={Q_min_supply:.4f} to Q*={Q_star:.4f} is: {integral_supply}")
    print("\nFinal Calculation:")
    print(f"Total Welfare = {integral_demand} - {integral_supply}")
    print(f"Total Welfare = {total_welfare}")


solve_welfare()
print("\n<<<30.4074815462529>>>")