import math
import numpy as np

def solve_hypercube_meeting_time():
    """
    This function calculates the expected time and variance for Alice and Bob to meet on a hypercube.
    The method involves analyzing the random walk of the relative position of Alice and Bob, Z_t = A_t - B_t (mod 2).
    The problem reduces to finding the hitting time to the origin for this new random walk. The state of this walk
    can be simplified to its Hamming distance from the origin. For an even dimension d, this distance is always even.

    Let e_k be the expected time to reach state 0 from distance k. We solve the system of linear equations
    (I - Q)e = 1, where Q is the transition probability matrix for the transient states {2, 4, ..., d}.
    
    Similarly, for the second moment f_k = E[T^2|start at k], we solve (I - Q)f = 1 + 2Qe.
    The variance from state k is f_k - (e_k)^2.
    """
    d = 14
    m = d // 2  # Number of transient states {k=2, 4, ..., 14}

    # Transition matrix Q for transient states {k=2, 4, ..., 2m}
    # States are indexed 0 to m-1, corresponding to k = 2, 4, ..., 14
    Q = np.zeros((m, m))
    
    for i in range(m):
        k = 2 * (i + 1)
        
        # Transition to k-2 (state i-1)
        if k > 0:
            prob_down = k * (k - 1) / (d**2)
            if i > 0:
                Q[i, i - 1] = prob_down

        # Transition to k+2 (state i+1)
        if k < d:
            prob_up = (d - k) * (d - k - 1) / (d**2)
            if i < m - 1:
                Q[i, i + 1] = prob_up
        else:
            prob_up = 0
            
        # Transition to k (state i) (staying in the same distance state)
        prob_stay_lazy = 1.0/d
        prob_stay_active = 2.0*k*(d-k) / d**2
        prob_stay = prob_stay_lazy + prob_stay_active
        Q[i, i] = prob_stay

    I = np.identity(m)
    A = I - Q
    
    # Solve for expected times e_k
    ones_vector = np.ones(m)
    e_vec = np.linalg.solve(A, ones_vector)
    
    # Expected time starting from distance d=14 is the last element
    e_d_14 = e_vec[-1]

    # Solve for second moments f_k
    b_f_vector = ones_vector + 2 * Q @ e_vec
    f_vec = np.linalg.solve(A, b_f_vector)
    
    # Second moment from distance d=14
    f_d_14 = f_vec[-1]

    var_d_14 = f_d_14 - e_d_14**2

    # E[X_15]
    # For d=15 (odd), Alice and Bob start on vertices of different parity
    # and always move to vertices of different parity. They can never meet.
    e_d_15_str = "inf"

    # Is it true that E[X_d] <= d/2 * d^d/d! for even d?
    def get_ex_d(dim):
        if dim % 2 != 0: return float('inf')
        m_d = dim // 2
        Q_d = np.zeros((m_d, m_d))
        for i in range(m_d):
            k = 2 * (i + 1)
            if i > 0:
                Q_d[i, i - 1] = k * (k - 1) / (dim**2)
            if i < m_d - 1:
                Q_d[i, i + 1] = (dim - k) * (dim - k - 1) / (dim**2)
            Q_d[i, i] = 1/dim + 2*k*(dim-k)/(dim**2)
        A_d = np.identity(m_d) - Q_d
        e_vec_d = np.linalg.solve(A_d, np.ones(m_d))
        return e_vec_d[-1]

    is_inequality_true = "yes"
    for d_test in [2, 4, 6, 8, 10]:
        try:
            lhs = get_ex_d(d_test)
            rhs = (d_test / 2.0) * (d_test**d_test / math.factorial(d_test))
            if lhs > rhs + 1e-9: # Add tolerance for float comparison
                is_inequality_true = "no"
                break
        except Exception:
            is_inequality_true = "Could not determine"
            break
            
    print(f"E[X_14]: {math.floor(e_d_14)}")
    print(f"D^2[X_14]: {math.floor(var_d_14)}")
    print(f"E[X_15]: {e_d_15_str}")
    print(f"Is E[X_d] <= d/2 * d^d/d! for even d? {is_inequality_true}")

solve_hypercube_meeting_time()