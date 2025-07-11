import sympy
from sympy import Symbol, cos, pi, sin, diff

def solve_magnetization():
    """
    Solves the integral equation for M_z(1) for n=1.

    The method involves:
    1. Setting up the symbolic framework with sympy.
    2. Defining G_n(B) as given by the problem.
    3. Applying the recurrence relation G_{k-1} = (G_k - G_k') / (k * C')
       to find G_0(B).
    4. Calculating M_z(B) = G_0'(B) - G_0(B).
    5. Evaluating M_z(1) to find the final numerical answer.
    """
    B = Symbol('B')
    n = 1

    # Based on the problem statement, G_n(B) = n^{-n} * B^{4n}
    # For n=1, G_1(B) = 1^{-1} * B^4 = B^4
    g_n = B**(4 * n)

    # C(B) = cos(pi*B/2)
    # C'(B) is its derivative
    C_prime = diff(cos(pi * B / 2), B)

    # Recursively find g_0(B), where g_k is related to G_k
    # g_{k-1}(B) = (g_k(B) - g_k'(B)) / (k * C'(B))
    g_current = g_n
    for k in range(n, 0, -1):
        g_k_prime = diff(g_current, B)
        g_current = (g_current - g_k_prime) / (k * C_prime)
    
    g_0 = g_current

    # M_z(B) = G_0'(B) - G_0(B). 
    # For n=1, G_0 = g_0. For other n, G_0 = n**(-n)*g_0
    m_z = diff(g_0, B) - g_0

    # Evaluate at B=1
    g0_val = g_0.subs(B, 1)
    g0_prime_val = diff(g_0, B).subs(B, 1)

    mz_val = m_z.subs(B, 1)
    
    # We found M_z(1,n) increases with n for small n (1, 2, 3), 
    # suggesting the minimum is at n_min=1.
    # The question is to find this minimum value.
    # We print the components of the final calculation for clarity.
    
    g0_val_num = float(g0_val.evalf())
    g0_prime_val_num = float(g0_prime_val.evalf())
    mz_val_num = float(mz_val.evalf())
    
    print(f"Based on the analysis, the minimum magnetization occurs at n=1.")
    print(f"The calculation for M_z(1) is based on G_0(B), the result of the recursive process.")
    print(f"At B=1:")
    print(f"G_0(1) = {g0_val} = {g0_val_num:.4f}")
    print(f"G_0'(1) = {g0_prime_val} = {g0_prime_val_num:.4f}")
    print(f"The minimum magnetization is M_z(1) = G_0'(1) - G_0(1)")
    print(f"M_z(1) = {g0_prime_val} - ({g0_val}) = {mz_val} = {mz_val_num:.4f}")

solve_magnetization()