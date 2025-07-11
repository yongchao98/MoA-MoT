import math
from scipy.optimize import root_scalar
from scipy.special import comb

def solve_equilibrium_p():
    """
    Solves for the equilibrium probability p in the described game.
    """
    
    # Parameters from the problem
    N = 8.0
    M = 3.0 * N
    M1 = M - 1.0

    def h(p):
        """
        Calculates the difference between the payoff of the discrete strategy
        and the payoff of the 2-race deviation strategy.
        The root of this function gives the equilibrium probability p.
        h(p) = U_discrete(p) - U_dev(p)
        """
        if p <= 0 or p > 1:
            return float('nan')
        
        # Payoff for the discrete strategy
        # U_discrete = (1 - (1 - p/N)**M) / (3 * p)
        # We use a check for small p to handle the limit correctly.
        term_d = (1 - p/N)
        if abs(p) < 1e-9:
             # Limit of U_discrete as p -> 0 is 1, using L'Hopital's Rule.
            U_discrete = 1.0
        else:
            U_discrete = (1 - term_d**M) / (3 * p)

        # Payoff for the deviation strategy (spreading over 2 races)
        term_dev1 = (1 - p/N)
        term_dev2 = (1 - 2*p/N)
        U_dev = 2 * (term_dev1)**(M1) - (term_dev2)**(M1)
        
        return U_discrete - U_dev

    # Find the root of the function h(p) in the interval (0, 1)
    # Based on analysis, h(p) is negative for small p and positive for p=1,
    # so a root must exist in (0,1).
    try:
        sol = root_scalar(h, bracket=[1e-6, 1.0])
        p = sol.root
    except ValueError:
        print("Error: Could not find a root in the given bracket. The model or analysis might be incorrect.")
        return

    # Print out the detailed steps as requested
    print("The equilibrium condition is found when the payoff from the discrete strategy equals the payoff from an alternative (deviating) strategy.")
    print(f"\nGame parameters: N = {int(N)} races, {int(M)} players.")
    print("Equation for equilibrium: U_discrete(p) = U_dev(p)\n")

    print("U_discrete(p) = (1 - (1 - p/8)^24) / (3*p)")
    print("U_dev(p)      = 2*(1 - p/8)^23 - (1 - 2*p/8)^23\n")

    print(f"Solving this equation numerically yields p = {p:.6f}\n")

    # Show that the payoffs are equal at this p
    U_d_val = (1 - (1 - p/N)**M) / (3 * p)
    U_dev_val = 2 * (1 - p/N)**(M1) - (1 - 2*p/N)**(M1)
    
    print("Verification at p = {:.6f}:".format(p))
    print("LHS = (1 - (1 - {:.6f}/8)^24) / (3*{:.6f}) = {:.6f}".format(p, p, U_d_val))
    print("RHS = 2*(1 - {:.6f}/8)^23 - (1 - 2*{:.6f}/8)^23 = {:.6f}".format(p, p, U_dev_val))
    print("\nThe payoffs are equal, confirming the equilibrium.\n")
    
    # Calculate the final requested value
    one_minus_p = 1 - p
    value_to_floor = 10000 * one_minus_p
    result = math.floor(value_to_floor)

    print("Final Calculation:")
    print(f"The value p is {p:.6f}")
    print(f"1 - p = {one_minus_p:.6f}")
    print(f"10000 * (1 - p) = {value_to_floor:.4f}")
    print(f"The floor of this value is {result}")
    
    # Final answer format
    print("\n<<<{}>>>".format(result))

if __name__ == '__main__':
    solve_equilibrium_p()