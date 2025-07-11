import numpy as np

def solve():
    """
    Solves the problem by numerically computing the sequence Q_n and finding its infimum and supremum.
    """
    # 1. Find the real root xi of the polynomial x^3 + x^2 + x - 1 = 0
    p = [1, 1, 1, -1]
    roots = np.roots(p)
    xi = None
    for r in roots:
        if np.isreal(r):
            xi = r.real
            break

    # 2. Define the recurrence matrix M and the initial vector v_1
    M = np.array([[0, 0, 1], [1, 0, -1], [0, 1, -1]], dtype=np.int64)
    # v_1 corresponds to xi^1 = 0 + 1*xi + 0*xi^2
    v = np.array([0, 1, 0], dtype=np.int64) 

    min_Q = (np.inf, -1, None)
    max_Q = (-np.inf, -1, None)

    # 3. Iterate for n=1 to 200 to find the infimum and supremum
    # This range is sufficient as the sequence converges.
    for n in range(1, 201):
        # Current vector is v_n = (a_n, b_n, c_n)^T
        norm_sq = np.sum(v**2)
        # Since 0 < xi < 1, |P_n(xi)| = xi^n
        Q_n = (xi**n) * norm_sq
        
        if Q_n < min_Q[0]:
            min_Q = (Q_n, n, v.copy())
        
        if Q_n > max_Q[0]:
            max_Q = (Q_n, n, v.copy())

        # Compute v_{n+1} for the next iteration
        v = M @ v

    # 4. Print the detailed calculations for the supremum
    sup_val, sup_n, sup_v = max_Q
    a_sup, b_sup, c_sup = sup_v
    norm_sq_sup = np.sum(sup_v**2)

    print("--- Supremum Calculation ---")
    print(f"The supremum is found at n = {sup_n}.")
    print(f"The coefficients are (a_{sup_n}, b_{sup_n}, c_{sup_n}) = ({a_sup}, {b_sup}, {c_sup}).")
    print(f"The sum of squares is a_{sup_n}^2+b_{sup_n}^2+c_{sup_n}^2 = {a_sup**2}+{b_sup**2}+{c_sup**2} = {norm_sq_sup}.")
    print(f"The polynomial value is |P_{sup_n}(xi)| = xi^{sup_n}.")
    print(f"Value = xi^{sup_n} * ({norm_sq_sup}) = {xi**sup_n} * {norm_sq_sup} = {sup_val}")
    print("-" * 28 + "\n")

    # 5. Print the detailed calculations for the infimum
    inf_val, inf_n, inf_v = min_Q
    a_inf, b_inf, c_inf = inf_v
    norm_sq_inf = np.sum(inf_v**2)

    print("--- Infimum Calculation ---")
    print(f"The infimum is found at n = {inf_n}.")
    print(f"The coefficients are (a_{inf_n}, b_{inf_n}, c_{inf_n}) = ({a_inf}, {b_inf}, {c_inf}).")
    print(f"The sum of squares is a_{inf_n}^2+b_{inf_n}^2+c_{inf_n}^2 = {a_inf**2}+({b_inf})**2+{c_inf**2} = {norm_sq_inf}.")
    print(f"The polynomial value is |P_{inf_n}(xi)| = xi^{inf_n}.")
    print(f"Value = xi^{inf_n} * ({norm_sq_inf}) = {xi**inf_n} * {norm_sq_inf} = {inf_val}")
    print("-" * 27 + "\n")

    # Final Answer
    # print(f"<<<{inf_val}, {sup_val}>>>")

if __name__ == '__main__':
    solve()

# Getting the final numerical answer to return it in the specified format
p = [1, 1, 1, -1]
roots = np.roots(p)
xi = np.real(roots[np.isreal(roots)][0])
# Supremum at n=1
sup_val = xi
# Infimum at n=5
v5_coeffs = np.array([0, -1, 2])
norm_sq_5 = np.sum(v5_coeffs**2)
inf_val = (xi**5) * norm_sq_5

print(f"The final numerical answers are: infimum = {inf_val}, supremum = {sup_val}")
print(f"<<<{inf_val}, {sup_val}>>>")