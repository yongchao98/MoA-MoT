import numpy as np

def calculate_lower_bound():
    """
    Calculates and prints the lower bound on the minimax risk R^*_n
    for a specific Gaussian location model instance.

    This function follows the plan outlined:
    1. Define a specific scenario (Gaussian location model with chosen parameters).
    2. Calculate delta, the separation parameter.
    3. Calculate the chi-squared divergence term.
    4. Use chi-squared to bound the total variation distance.
    5. Compute the final lower bound on the minimax risk.
    """

    # --- Step 1: Define a concrete scenario ---
    # We instantiate the general problem with a 1D Gaussian location model.
    # The distributions are P_j = N(mu_j, 1) and P_0 = N(0, 1).
    # The parameter of interest is the mean, theta_j = mu_j.
    # We choose specific values for the problem parameters.
    
    # Number of i.i.d. samples
    n = 10
    
    # Number of alternative hypotheses
    N = 4
    
    # Means for the N alternative hypotheses
    mus = np.array([0.1, -0.1, 0.2, -0.2])
    
    # The metric rho is the Euclidean distance |a - b|.
    # The function Phi is chosen to be the identity, Phi(x) = x.
    def Phi(x):
        return x

    print("--- Scenario Parameters ---")
    print(f"Number of samples (n): {n}")
    print(f"Number of alternative hypotheses (N): {N}")
    print(f"Means of alternative distributions (mu_j): {mus}")
    print("Base distribution P_0 is N(0, 1)")
    print("Metric rho is Euclidean distance, Phi is identity function.")
    print("-" * 30 + "\n")

    # --- Step 2: Calculate delta ---
    # delta = min_{j in {1,...,N}} rho(theta(P_0), theta(P_j))
    # Here, theta(P_j) = mu_j and theta(P_0) = 0. rho is |.|.
    delta = np.min(np.abs(mus))
    
    # The term Phi(delta / 2) from the bound formula
    phi_delta_half = Phi(delta / 2.0)

    print("--- Intermediate Calculation Steps ---")
    print("Equation for the lower bound: R*_n >= (Phi(delta/2) / 2) * (1 - TV_distance_bound)")
    print(f"1. Delta (min distance to null): {delta}")
    print(f"2. Phi(delta / 2): {phi_delta_half}")

    # --- Step 3: Calculate the Chi-squared divergence term ---
    # For N(mu, 1) vs N(0, 1), the term E_P0[(dPj/dP0)(dPk/dP0)] is exp(mu_j * mu_k)
    # chi_squared = (1/N^2) * sum_{j,k} (E_P0[...])^n - 1
    #             = (1/N^2) * sum_{j,k} exp(n * mu_j * mu_k) - 1
    sum_term = 0
    for mu_j in mus:
        for mu_k in mus:
            sum_term += np.exp(n * mu_j * mu_k)
    
    chi_squared = (sum_term / (N**2)) - 1
    
    print(f"3. Chi-squared divergence term: {chi_squared:.5f}")

    # --- Step 4: Bound the Total Variation (TV) distance ---
    # We use the inequality: d_TV(P, Q) <= 0.5 * sqrt(chi_squared(Q, P))
    # We use max(0, chi_squared) to avoid domain errors if chi_squared is slightly negative due to precision.
    tv_distance_bound = 0.5 * np.sqrt(max(0, chi_squared))

    print(f"4. Upper bound on the Total Variation distance: {tv_distance_bound:.5f}")

    # --- Step 5: Compute the final lower bound on R^*_n ---
    # R^*_n >= (Phi(delta/2) / 2) * (1 - TV_distance_bound)
    lower_bound = (phi_delta_half / 2.0) * (1 - tv_distance_bound)
    # The bound is trivial (<= 0) if the hypotheses are easily distinguishable.
    lower_bound = max(0, lower_bound)
    
    print("-" * 30 + "\n")
    print("--- Final Result ---")
    print(f"The final equation is: {phi_delta_half / 2.0:.4f} * (1 - {tv_distance_bound:.4f})")
    print(f"Calculated tightest lower bound on R^*_n: {lower_bound:.5f}")
    
    return lower_bound

if __name__ == '__main__':
    result = calculate_lower_bound()
    # The final answer format is specified in the prompt.
    # print(f"\n<<<{result:.5f}>>>")