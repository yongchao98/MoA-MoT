import numpy as np

def solve_polynomial_inf_sup():
    """
    Calculates the infimum and supremum of the given expression.
    """
    # Find the real root xi of the polynomial x^3 + x^2 + x - 1 = 0
    roots = np.roots([1, 1, 1, -1])
    xi = roots[np.isreal(roots)][0].real

    # Matrix for the recurrence relation of coefficients (a_n, b_n, c_n)
    M = np.array([[0, 0, 1], [1, 0, -1], [0, 1, -1]])

    # Initial coefficient vector for n=1, P_1(x) = x
    v = np.array([0, 1, 0])

    q_values = []
    v_values = []

    # We compute for a range of n to find the min/max values.
    # The sequence converges to a limit cycle, so the extrema will be found for small n.
    # A range up to n=50 is sufficient to observe this behavior.
    for n in range(1, 51):
        # Store the coefficient vector for the current n
        v_values.append(v)
        
        # Calculate Q_n = xi^n * (a_n^2 + b_n^2 + c_n^2)
        q_n = (xi**n) * np.sum(v**2)
        q_values.append(q_n)
        
        # Update the coefficient vector for n+1
        v = M @ v

    # Find the infimum and supremum from the computed list
    inf_val = min(q_values)
    sup_val = max(q_values)
    
    n_inf = q_values.index(inf_val) + 1
    n_sup = q_values.index(sup_val) + 1
    
    v_inf = v_values[n_inf - 1]
    v_sup = v_values[n_sup - 1]

    print("Infimum calculation:")
    print(f"The infimum is found at n = {n_inf}")
    print(f"Coefficients (a_n, b_n, c_n) are ({v_inf[0]}, {v_inf[1]}, {v_inf[2]})")
    print(f"The value is |P_{n_inf}(xi)| * (a_{n_inf}^2 + b_{n_inf}^2 + c_{n_inf}^2) = xi^{n_inf} * ({v_inf[0]**2} + {v_inf[1]**2} + {v_inf[2]**2})")
    print(f"inf_n |P_n(xi)| (a_n^2+b_n^2+c_n^2) = {inf_val}")
    
    print("\n" + "="*30 + "\n")

    print("Supremum calculation:")
    print(f"The supremum is found at n = {n_sup}")
    print(f"Coefficients (a_n, b_n, c_n) are ({v_sup[0]}, {v_sup[1]}, {v_sup[2]})")
    print(f"The value is |P_{n_sup}(xi)| * (a_{n_sup}^2 + b_{n_sup}^2 + c_{n_sup}^2) = xi^{n_sup} * ({v_sup[0]**2} + {v_sup[1]**2} + {v_sup[2]**2})")
    print(f"sup_n |P_n(xi)| (a_n^2+b_n^2+c_n^2) = {sup_val}")

solve_polynomial_inf_sup()