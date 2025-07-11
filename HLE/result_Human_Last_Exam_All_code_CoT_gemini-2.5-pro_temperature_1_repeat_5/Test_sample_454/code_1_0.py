import numpy as np
from scipy.optimize import fsolve

def solve_adsorption_system():
    """
    Solves the multilayer adsorption problem using mean-field theory.
    """
    # Given parameters in dimensionless form
    beta_mu = 0.1
    beta_epsilon = -1 / (2 * np.pi)

    # --- Step 1: Solve for the conditional occupancy of upper layers <n_L> ---
    # The equation is: ln(<n_L> / (1 - <n_L>)) = beta*(mu_L - E_int,L)
    # With mu_L = 0 and E_int,L = (8*<n_L> + 4)*epsilon, this becomes:
    # ln(<n_L> / (1 - <n_L>)) = -(8*<n_L> + 4)*beta_epsilon
    
    def equation_nL(nL):
        if nL <= 0 or nL >= 1:
            return 1e9 # Return a large number for invalid inputs
        return np.log(nL / (1 - nL)) + (8 * nL + 4) * beta_epsilon

    # Initial guess for the solution
    initial_guess_nL = 0.8
    nL_solution = fsolve(equation_nL, initial_guess_nL)[0]
    
    print(f"The conditional occupancy of upper layers is <n_L> = {nL_solution:.3f}")

    # --- Step 2: Solve for the conditional occupancy of the first layer <n_1> ---
    # The equation is: ln(<n_1> / (1 - <n_1>)) = beta*mu - beta*E_int,1
    # E_int,1 = (4*<n_1> + 4*<n_L>)*epsilon
    # ln(<n_1> / (1 - <n_1>)) = beta_mu - (4*<n_1> + 4*<n_L>)*beta_epsilon
    
    def equation_n1(n1, nL):
        if n1 <= 0 or n1 >= 1:
            return 1e9 # Return a large number for invalid inputs
        return np.log(n1 / (1 - n1)) - (beta_mu - (4 * n1 + 4 * nL) * beta_epsilon)

    # Initial guess for the solution
    initial_guess_n1 = 0.7
    # Use a lambda function to pass the solved nL_solution as a fixed parameter
    n1_solution = fsolve(lambda n1: equation_n1(n1, nL_solution), initial_guess_n1)[0]
    
    print(f"The conditional occupancy of the first layer is <n_1> = {n1_solution:.3f}")

    # --- Step 3: Calculate the total average occupancy <n> ---
    # <n> = <n_1> / (1 - <n_L>)
    
    n_total = n1_solution / (1 - nL_solution)
    
    print("\nThe total average occupancy is calculated using the formula: <n> = <n_1> / (1 - <n_L>)")
    print(f"{n_total:.3f} = {n1_solution:.3f} / (1 - {nL_solution:.3f})")
    
    print(f"\nThe final average occupancy per site <n> is: {n_total:.3f}")
    
    return n_total

if __name__ == '__main__':
    final_answer = solve_adsorption_system()
    print(f"\n<<<{final_answer:.3f}>>>")
