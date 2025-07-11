import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def solve_and_print_welfare():
    """
    This function solves for market equilibrium, calculates consumer surplus,
    producer surplus, and total welfare, and prints the results.
    """
    # 1. Define the supply and demand functions
    def demand(Q):
        return 18 * np.exp(-np.arctan(Q))

    def supply(Q):
        # The function is defined for Q^3 - 2 > 0
        if Q**3 - 2 <= 0:
            return np.inf  # Return infinity if outside the domain
        return np.log(Q**3 - 2)

    # 2. Find the market equilibrium
    # We need to solve Supply(Q) = Demand(Q), or Supply(Q) - Demand(Q) = 0
    def equilibrium_equation(Q):
        # fsolve expects a function that returns 0 at the root.
        # Handle the domain of the supply function.
        if Q**3 - 2 <= 0:
            return 1e10  # Return a large number if outside the domain
        return supply(Q) - demand(Q)

    # Initial guess for Q. From analysis, the root is between 4 and 5.
    initial_guess_q = 4.5
    # Use fsolve to find the root Q_E
    q_equilibrium = fsolve(equilibrium_equation, initial_guess_q)[0]
    # Calculate P_E using the demand function
    p_equilibrium = demand(q_equilibrium)

    # 3. Calculate Consumer and Producer Surplus using numerical integration
    
    # Consumer Surplus calculation
    # CS = Integral(Demand) from 0 to Q_E  -  (P_E * Q_E)
    integral_demand, _ = quad(demand, 0, q_equilibrium)
    consumer_surplus = integral_demand - p_equilibrium * q_equilibrium

    # Producer Surplus calculation
    # PS = (P_E * Q_E) - Integral(Supply) from Q_min to Q_E
    # The supply curve starts where Q^3 - 2 > 0, so Q > 2^(1/3)
    q_supply_starts = 2**(1/3)
    # Define a simple supply function for integration without the conditional check
    def supply_for_quad(Q):
        return np.log(Q**3 - 2)
    integral_supply, _ = quad(supply_for_quad, q_supply_starts, q_equilibrium)
    producer_surplus = p_equilibrium * q_equilibrium - integral_supply

    # 4. Calculate Total Welfare
    total_welfare = consumer_surplus + producer_surplus

    # 5. Print the results in a clear, step-by-step format
    print("--- Market Equilibrium ---")
    print(f"The equilibrium equation is: ln(Q\u00b3 - 2) = 18 * e^(-arctan(Q))")
    print(f"Solving this numerically, we find the equilibrium point:")
    print(f"Equilibrium Quantity (Q_E) = {q_equilibrium:.4f}")
    print(f"Equilibrium Price (P_E) = {p_equilibrium:.4f}\n")

    print("--- Welfare Calculation ---")
    # Print Consumer Surplus equation and result
    print("Consumer Surplus (CS) = \u222b[0 to Q_E]D(Q)dQ - (P_E * Q_E)")
    print(f"CS = {integral_demand:.4f} - ({p_equilibrium:.4f} * {q_equilibrium:.4f})")
    print(f"CS = {consumer_surplus:.4f}\n")

    # Print Producer Surplus equation and result
    print("Producer Surplus (PS) = (P_E * Q_E) - \u222b[2\u207d\u00b9\u2044\u00b3\u207e to Q_E]S(Q)dQ")
    print(f"PS = ({p_equilibrium:.4f} * {q_equilibrium:.4f}) - {integral_supply:.4f}")
    print(f"PS = {producer_surplus:.4f}\n")
    
    # Print Total Welfare equation and final result
    print("--- Total Welfare ---")
    print("Total Welfare (TW) = Consumer Surplus + Producer Surplus")
    print(f"TW = {consumer_surplus:.4f} + {producer_surplus:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")
    
    return total_welfare

# Execute the function and capture the final answer
final_answer = solve_and_print_welfare()
print(f"\n<<< {final_answer:.4f} >>>")
