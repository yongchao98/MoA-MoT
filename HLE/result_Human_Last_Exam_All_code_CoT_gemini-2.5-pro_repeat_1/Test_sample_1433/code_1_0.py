import numpy as np
import math

def solve_hypercube_meeting_time(d):
    """
    Calculates the expected time and variance for meeting on a d-hypercube for two
    people starting at opposite corners.
    d must be an even integer >= 2.
    """
    if d % 2 != 0:
        raise ValueError("d must be an even integer for this function.")
    if d < 2:
        raise ValueError("d must be >= 2.")

    if d == 2:
        # For d=2, the meeting probability at each step is 1/2.
        # This is a geometric distribution with p=1/2.
        # E = 1/p = 2. Var = (1-p)/p^2 = (1/2)/(1/4) = 2.
        return 2.0, 2.0
    
    # The states of the Markov chain on distance k are {0, 2, ..., d}.
    # The transient states are {2, 4, ..., d-2}.
    num_transient_states = d // 2 - 1
    
    # A is the submatrix of (I-P) for transient states.
    # We solve A*x = b.
    A = np.zeros((num_transient_states, num_transient_states))
    b_E = np.ones(num_transient_states)
    
    # Populate the matrix A for the system on E_k
    for i in range(num_transient_states):
        k = 2 * (i + 1)
        p_km2 = (k * (k - 1)) / (d * d)
        p_kp2 = ((d - k) * (d - k - 1)) / (d * d)
        p_k = 1.0 - p_km2 - p_kp2 # Probability of staying at distance k
        
        # From E_k = 1 + p_k*E_k + p_km2*E_{k-2} + p_kp2*E_{k+2}
        # (1-p_k)E_k - p_km2*E_{k-2} - p_kp2*E_{k+2} = 1
        A[i, i] = 1.0 - p_k
        if i > 0:
            A[i, i-1] = -p_km2
        if i < num_transient_states - 1:
            A[i, i+1] = -p_kp2

    # The equation for k=d-2 depends on E_d.
    # We use the relation E_d = E_{d-2} + d/(d-1).
    k_last = d - 2
    p_last_d = ((d - k_last) * (d - k_last - 1)) / (d * d)
    # The term -p_last_d * E_d becomes -p_last_d * (E_{d-2} + d/(d-1))
    A[num_transient_states-1, num_transient_states-1] -= p_last_d
    b_E[num_transient_states-1] += p_last_d * d / (d - 1)

    # Solve for E_k for k in {2, ..., d-2}
    E_vec = np.linalg.solve(A, b_E)
    E_dm2 = E_vec[-1]
    E_d_val = E_dm2 + d / (d - 1)

    # --- Variance Calculation ---
    # Using G = (2N-I)E for the transient states, where N = inv(A) and G is the vector of 2nd moments.
    N_matrix = np.linalg.inv(A)
    G_vec = (2 * N_matrix - np.identity(num_transient_states)) @ E_vec
    G_dm2 = G_vec[-1]

    # We need G_d. The relation is derived from G_d = 1 + P_dd(G_d+2E_d+1) + P_d,d-2(G_{d-2}+2E_{d-2}+1)
    P_d_dm2 = (d * (d - 1)) / (d * d)
    P_dd = d / (d * d)
    
    # (1-P_dd)G_d - P_d,d-2*G_{d-2} = 1 + P_dd(2E_d+1) + P_d,d-2(2E_{d-2}+1)
    # P_d,d-2 * (G_d - G_{d-2}) = 2 + 2*P_dd*E_d + 2*P_d,d-2*E_{d-2}
    rhs_G = 2 + 2 * P_dd * E_d_val + 2 * P_d_dm2 * E_dm2
    G_d_val = G_dm2 + rhs_G / P_d_dm2
    
    V_d_val = G_d_val - E_d_val**2
    
    return E_d_val, V_d_val

# --- Main Calculations ---

# Problem 1: E[X_14] and D^2[X_14]
E14, V14 = solve_hypercube_meeting_time(14)
print(f"The expected time E[X_14] is {E14:.4f}, its integer part is: {math.floor(E14)}")
print(f"The variance D^2[X_14] is {V14:.4f}, its integer part is: {math.floor(V14)}")

# Problem 2: E[X_15]
E15_str = "inf"
print(f"The expected time E[X_15] is: {E15_str}")

# Problem 3: Inequality check
inequality_answer = "yes"
print(f"Is it true that for even d, E[X_d] <= (d/2) * d^d / d!? {inequality_answer}")

# Final Answer Block
print(f"<<<{math.floor(E14)}, {math.floor(V14)}, {E15_str}, {inequality_answer}>>>")