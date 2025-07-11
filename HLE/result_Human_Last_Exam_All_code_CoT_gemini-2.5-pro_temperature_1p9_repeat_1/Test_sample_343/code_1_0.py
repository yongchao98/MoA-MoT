import numpy as np

def solve_for_clique_count(k, T_values, P_values):
    """
    Solves for the number of k-cliques given the results from PCount.
    The equation is P(T) = c1 * T^k + c2 * T^(2k)

    Args:
        k (int): The parameter k.
        T_values (list): A list of T values used, e.g., [1, 2].
        P_values (list): The corresponding results from PCount(G_T, k).
    """
    # We need to solve a system of linear equations for the coefficients.
    # For k=2, we have:
    # P(T1) = c1*T1^2 + c2*T1^4
    # P(T2) = c1*T2^2 + c2*T2^4
    #
    # In matrix form Ax = b:
    # [[T1^k, T1^(2k)], [T2^k, T2^(2k)]] [c1, c2]' = [P1, P2]'

    if len(T_values) < 2 or len(P_values) < 2:
        print("Need at least 2 points to solve for 2 coefficients.")
        return

    T1, T2 = T_values[0], T_values[1]
    P1, P2 = P_values[0], P_values[1]

    # Constructing the matrix A
    A = np.array([
        [T1**k, T1**(2*k)],
        [T2**k, T2**(2*k)]
    ])

    # Constructing the vector b
    b = np.array([P1, P2])

    try:
        # Solve the system Ax = b for x = [c1, c2]
        coeffs = np.linalg.solve(A, b)
        num_cliques = coeffs[0]
        
        # We model the equation P(T) = N_clique * T^k + (N_im + N_ib) * T^2k
        # Outputting the numbers from the equation for one of the data points, T1.
        print(f"For T = {T1}, the equation is: {P1} = {coeffs[0]:.1f} * {T1}^{k} + {coeffs[1]:.1f} * {T1}^{2*k}")
        print("Each number in the final equation:")
        print(f"Total Count: {P1}")
        print(f"Coefficient 1 (Number of cliques): {coeffs[0]:.1f}")
        print(f"Exponent 1: {k}")
        print(f"Coefficient 2 (Sum of other counts): {coeffs[1]:.1f}")
        print(f"Exponent 2: {2*k}")

    except np.linalg.LinAlgError:
        print("Could not solve the system. The matrix A may be singular.")


# --- Example Usage ---
# Suppose for a graph G and k=2, we construct G_T for T=1 and T=2.
# We run PCount and get (hypothetical) results:
# For T=1, PCount(G_1, 2) = 15  (i.e., N_clique + N_im + N_ib = 15)
# For T=2, PCount(G_2, 2) = 116 (i.e., N_clique*4 + (N_im + N_ib)*16 = 116)
k_param = 2
T_vals = [1, 2]
P_vals = [15, 116]

solve_for_clique_count(k_param, T_vals, P_vals)
