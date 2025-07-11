import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Step 1: Define the supply and demand functions
def supply_func(Q):
    """Defines the supply function P = ln(Q^3 - 2)"""
    # The function is only defined for Q^3 - 2 > 0
    if Q**3 - 2 <= 0:
        return np.inf  # Return a very large number for invalid inputs to guide the solver
    return np.log(Q**3 - 2)

def demand_func(Q):
    """Defines the demand function P = 18 * e^(-arctan(Q))"""
    return 18 * np.exp(-np.arctan(Q))

# Step 2: Find the market equilibrium
def find_equilibrium():
    """Finds the equilibrium quantity (Q_E) and price (P_E)"""
    # Define an equation representing the difference between demand and supply
    def equilibrium_equation(Q):
        return demand_func(Q) - supply_func(Q)

    # Provide an initial guess for the solver. A quick plot or trial and error suggests Q is around 4.5.
    initial_guess_Q = 4.5
    # Use fsolve to find the root of the equilibrium equation, which is Q_E
    Q_E_solution = fsolve(equilibrium_equation, initial_guess_Q)
    Q_E = Q_E_solution[0]
    
    # Calculate P_E using the demand function
    P_E = demand_func(Q_E)
    
    return Q_E, P_E

# Step 3: Calculate total welfare
def calculate_total_welfare(Q_E):
    """Calculates total welfare by integrating the area between the curves"""
    # The lower bound for the supply integral is where Q^3 - 2 = 0 => Q = 2^(1/3)
    supply_q_min = 2**(1/3)

    # Calculate the definite integral of the demand function from 0 to Q_E
    integral_demand, _ = quad(demand_func, 0, Q_E)

    # Calculate the definite integral of the supply function from its starting point to Q_E
    integral_supply, _ = quad(supply_func, supply_q_min, Q_E)

    # Total welfare is the difference between the two integrals
    total_welfare = integral_demand - integral_supply
    
    return integral_demand, integral_supply, total_welfare, supply_q_min

# --- Main execution ---
if __name__ == "__main__":
    # Find equilibrium
    Q_E, P_E = find_equilibrium()
    
    # Calculate welfare and its components
    integral_demand, integral_supply, total_welfare, supply_q_min = calculate_total_welfare(Q_E)
    
    # Step 4: Display the results
    print(f"First, we find the market equilibrium point (Q_E, P_E):")
    print(f"Equilibrium Quantity (Q_E): {Q_E:.4f}")
    print(f"Equilibrium Price (P_E): {P_E:.4f}\n")
    
    print(f"Next, we calculate the total welfare using integration.")
    print(f"Total Welfare = (Integral of Demand) - (Integral of Supply)")
    print(f"The equation with values is:")
    print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
    print(f"\nFinal Total Welfare = {total_welfare:.4f}")

    # Output the final answer in the specified format
    print(f"<<<${total_welfare:.4f}>>>")