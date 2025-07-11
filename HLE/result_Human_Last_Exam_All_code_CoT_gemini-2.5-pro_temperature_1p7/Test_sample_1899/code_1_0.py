import numpy as np

def solve_polynomial_sequence():
    """
    This function calculates the infimum and supremum of the sequence
    |P_n(xi)| * (a_n^2 + b_n^2 + c_n^2)
    as described in the problem.
    """
    
    # 1. Find the real root of the polynomial x^3 + x^2 + x - 1 = 0
    coeffs = [1, 1, 1, -1]
    roots = np.roots(coeffs)
    
    xi = 0.0
    for r in roots:
        if abs(r.imag) < 1e-12:
            xi = r.real
            break
            
    # 2. Set up the recurrence relation
    # Initial vector for n=1 (P_1(x) = x => a_1=0, b_1=1, c_1=0)
    v_n = np.array([0.0, 1.0, 0.0])

    # Recurrence matrix A
    A = np.array([
        [0, 0, 1],
        [1, 0, -1],
        [0, 1, -1]
    ], dtype=float)

    q_values = []
    # Maximum n for iteration. 50 steps is sufficient to observe the behavior.
    max_n = 50 
    
    # 3. Compute the sequence Q_n
    for n in range(1, max_n + 1):
        # Coefficients are given by the vector v_n = (a_n, b_n, c_n)
        a_n, b_n, c_n = v_n
        
        # Calculate norm squared: a_n^2 + b_n^2 + c_n^2
        norm_sq = np.dot(v_n, v_n)
        
        # P_n(xi) = xi^n
        p_n_xi_val = xi**n
        
        # Calculate Q_n
        q_n = abs(p_n_xi_val) * norm_sq
        q_values.append(q_n)
        
        # Update v_n for the next iteration: v_{n+1} = A @ v_n
        v_n = A @ v_n

    # 4. Determine infimum and supremum from the computed values
    inf_q = min(q_values)
    sup_q = max(q_values)
    
    inf_n = q_values.index(inf_q) + 1
    sup_n = q_values.index(sup_q) + 1

    print(f"The unique real root is xi ≈ {xi}")
    print("-" * 30)

    print(f"The supremum is found at n = {sup_n}.")
    v_sup = np.linalg.matrix_power(A, sup_n - 1) @ np.array([0.0, 1.0, 0.0])
    a_sup, b_sup, c_sup = v_sup
    print(f"For n={sup_n}: a_{sup_n}={a_sup:.0f}, b_{sup_n}={b_sup:.0f}, c_{sup_n}={c_sup:.0f}")
    print(f"sup |P_n(xi)|(a_n^2+b_n^2+c_n^2) = |{xi}^{sup_n}| * ({a_sup:.0f}^2 + {b_sup:.0f}^2 + {c_sup:.0f}^2) = {sup_q}")
    print("-" * 30)
    
    print(f"The infimum is found at n = {inf_n}.")
    v_inf = np.linalg.matrix_power(A, inf_n - 1) @ np.array([0.0, 1.0, 0.0])
    a_inf, b_inf, c_inf = v_inf
    print(f"For n={inf_n}: a_{inf_n}={a_inf:.0f}, b_{inf_n}={b_inf:.0f}, c_{inf_n}={c_inf:.0f}")
    print(f"inf |P_n(xi)|(a_n^2+b_n^2+c_n^2) = |{xi}^{inf_n}| * ({a_inf:.0f}^2 + {b_inf:.0f}^2 + {c_inf:.0f}^2) = {inf_q}")

solve_polynomial_sequence()