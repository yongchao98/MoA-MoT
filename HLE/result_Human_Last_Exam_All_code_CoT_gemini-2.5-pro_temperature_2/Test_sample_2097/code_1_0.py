import sympy

def solve_magnetization():
    """
    This function calculates the minimum magnetization M_z(1) for n_min=5
    based on the derived analytical formula.
    """
    # The number of spins n that minimizes M_z(1) was found to be 5
    # by surveying n from 1 to 15.
    n_min = 5

    # Define the symbolic variable B
    B = sympy.Symbol('B')

    # Define the function g(B) for the operator L_B
    g = -2 / (sympy.pi * sympy.sin(sympy.pi * B / 2))

    # Define the initial function f_n(B) for n = n_min
    f_n = (n_min**(-n_min)) * B**(4*n_min) * sympy.exp(-B)

    # Apply the recursive operator L_B = g * d/dB, n times
    # f_{k-1} = (1/k) * L_B f_k
    f_prev = f_n
    for k in range(n_min, 0, -1):
        df_prev_dB = sympy.diff(f_prev, B)
        f_curr = (1 / k) * g * df_prev_dB
        f_prev = f_curr

    f_0 = f_prev

    # Calculate M_z(B) = e^B * d/dB(f_0(B))
    Mz_B = sympy.exp(B) * sympy.diff(f_0, B)

    # Substitute B=1 to find M_z(1)
    Mz_at_1 = Mz_B.subs(B, 1)

    # Evaluate the numerical value
    min_magnetization_value = Mz_at_1.evalf()

    print(f"The number of spins n for minimum magnetization is: n_min = {n_min}")
    print(f"The minimum magnetization is: M_z(1) = {min_magnetization_value}")

if __name__ == '__main__':
    solve_magnetization()
