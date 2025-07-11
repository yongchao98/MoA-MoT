import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """
    
    # 1. Define supply and demand functions
    def demand(q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(q))

    def supply(q):
        """Supply function: P = ln(Q^3 - 2)"""
        # The function is only defined for q^3 - 2 > 0.
        # This check prevents math domain errors during numerical solving.
        # If q is an array, we handle it element-wise.
        q = np.asarray(q)
        # Create a result array of the same shape, initialized to infinity.
        result = np.full(q.shape, np.inf)
        # Identify valid q values
        valid_mask = q**3 > 2
        # Compute log for valid q values
        result[valid_mask] = np.log(q[valid_mask]**3 - 2)
        # If the input was a scalar, return a scalar
        return result.item() if q.ndim == 0 else result

    # 2. Find market equilibrium (Q_e, P_e)
    # We solve for the quantity Q where Demand(Q) = Supply(Q)
    def equilibrium_equation(q):
        return demand(q) - supply(q)

    # The supply curve is defined for Q > 2^(1/3) ~= 1.26.
    # We use an initial guess greater than this value for the solver.
    initial_guess_q = 2.0
    # Use fsolve to find the equilibrium quantity Q_e
    Q_e = fsolve(equilibrium_equation, initial_guess_q)[0]
    
    # Calculate the equilibrium price P_e using the demand function
    P_e = demand(Q_e)

    # 3. Calculate total welfare
    # Total Welfare (TW) = Consumer Surplus (CS) + Producer Surplus (PS)
    # TW = (Integral[0 to Qe] D(Q) dQ - Pe*Qe) + (Pe*Qe - Integral[Q_min to Qe] S(Q) dQ)
    # TW = Integral[0 to Qe] D(Q) dQ - Integral[Q_min to Qe] S(Q) dQ

    # Calculate the integral of the demand function from 0 to Q_e
    integral_demand, _ = quad(demand, 0, Q_e)
    
    # Define the lower bound for the supply curve integral
    q_min_supply = (2)**(1/3)
    
    # Calculate the integral of the supply function from its start point to Q_e
    integral_supply, _ = quad(supply, q_min_supply, Q_e)
    
    # Calculate the final total welfare
    total_welfare = integral_demand - integral_supply
    
    # 4. Output the results
    print(f"Finding the market equilibrium...")
    print(f"Equilibrium Quantity (Q_e): {Q_e:.4f}")
    print(f"Equilibrium Price (P_e): {P_e:.4f}\n")
    
    print("Calculating Total Welfare as the area between the demand and supply curves...")
    print("Total Welfare = ∫ D(Q) dQ - ∫ S(Q) dQ")
    
    # Print the final equation with calculated values
    final_equation_str_1 = f"Total Welfare = ∫[from 0 to {Q_e:.4f}] (18e^(-arctan(Q))) dQ - ∫[from {q_min_supply:.4f} to {Q_e:.4f}] (ln(Q³-2)) dQ"
    final_equation_str_2 = f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}"
    final_equation_str_3 = f"Total Welfare = {total_welfare:.4f}"
    
    print(final_equation_str_1)
    print(final_equation_str_2)
    print(final_equation_str_3)

if __name__ == "__main__":
    main()
<<<10.047318359419137>>>