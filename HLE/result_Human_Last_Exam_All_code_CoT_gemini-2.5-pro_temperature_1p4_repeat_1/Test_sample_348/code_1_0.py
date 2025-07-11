import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def solve_and_calculate_welfare():
    """
    This function defines the supply and demand curves, finds the market equilibrium,
    and calculates the consumer surplus, producer surplus, and total welfare.
    """
    
    # 1. Define the supply and demand functions
    # Demand: P = 18 * exp(-arctan(Q))
    def demand_func(q):
        return 18 * np.exp(-np.arctan(q))

    # Supply: P = ln(Q^3 - 2)
    def supply_func(q):
        # The function is defined for q^3 - 2 > 0.
        # Handle values outside the domain to assist the numerical solver.
        # We use a small positive number inside log to avoid -inf.
        arg = q**3 - 2
        if arg <= 0:
            return np.inf  # Return a large number to repel the solver from this region
        return np.log(arg)

    # 2. Find the market equilibrium
    # We solve for Q where Demand(Q) = Supply(Q), or Demand(Q) - Supply(Q) = 0
    def equilibrium_equation(q):
        return demand_func(q) - supply_func(q)

    # The supply curve is defined for Q > 2^(1/3) ~= 1.26.
    # An initial guess must be in this domain.
    initial_guess = 2.0
    
    # Use fsolve to find the root of the equilibrium equation
    q_equilibrium_array = fsolve(equilibrium_equation, initial_guess)
    q_equilibrium = q_equilibrium_array[0]

    # Calculate the equilibrium price using the demand function
    p_equilibrium = demand_func(q_equilibrium)

    # 3. Calculate Consumer and Producer Surplus
    
    # Consumer Surplus (CS) = integral from 0 to Q_E of (Demand(Q) - P_E) dQ
    cs_integrand = lambda q: demand_func(q) - p_equilibrium
    consumer_surplus, _ = quad(cs_integrand, 0, q_equilibrium)

    # The supply curve is valid for Q > 2^(1/3)
    q_min_supply = np.cbrt(2)
    
    # Producer Surplus (PS) = integral from Q_min to Q_E of (P_E - Supply(Q)) dQ
    ps_integrand = lambda q: p_equilibrium - supply_func(q)
    producer_surplus, _ = quad(ps_integrand, q_min_supply, q_equilibrium)
    
    # 4. Calculate Total Welfare
    total_welfare = consumer_surplus + producer_surplus

    # 5. Display the results
    print("--- Market Equilibrium ---")
    print(f"Equilibrium Quantity (Q_E): {q_equilibrium:.4f}")
    print(f"Equilibrium Price (P_E): {p_equilibrium:.4f}")
    print("\n--- Welfare Calculation ---")
    print(f"Consumer Surplus (CS): {consumer_surplus:.4f}")
    print(f"Producer Surplus (PS): {producer_surplus:.4f}")
    print("\n--- Total Welfare ---")
    print(f"Total Welfare = CS + PS")
    print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f}")
    print(f"Final Total Welfare = {total_welfare:.4f}")

# Execute the function
solve_and_calculate_welfare()