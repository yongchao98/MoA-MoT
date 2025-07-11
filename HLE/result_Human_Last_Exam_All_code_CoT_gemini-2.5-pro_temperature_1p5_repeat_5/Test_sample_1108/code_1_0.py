import sympy as sp

def solve_linearization():
    """
    This function calculates the coefficients of the linearized Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # 1. Define symbolic variables and parameters
    S, F = sp.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # 2. Define the system of differential equations
    # dS/dt = g(S, F)
    # dF/dt = k(S, F)
    g = S * (h - m * S / F)
    k = F * (a - b * F - c * S)

    # 3. Find the non-trivial equilibrium point (S_e, F_e)
    # We solve the system g=0, k=0 for S > 0, F > 0.
    # From g=0 (and S>0), we get h - m*S/F = 0  => S = h*F/m
    # Substitute S into k=0 (and F>0): a - b*F - c*(h*F/m) = 0
    # Solving for F: F_e = a / (b + c*h/m)
    F_e_val = a / (b + c * h / m)
    # Solving for S: S_e = h*F_e/m
    S_e_val = (h * F_e_val) / m
    
    equilibrium_point = {S: S_e_val, F: F_e_val}

    # 4. Compute the Jacobian matrix J = [[∂g/∂S, ∂g/∂F], [∂k/∂S, ∂k/∂F]]
    J = sp.Matrix([[sp.diff(g, S), sp.diff(g, F)],
                   [sp.diff(k, S), sp.diff(k, F)]])

    # 5. Evaluate the Jacobian at the equilibrium point to find the 'a' coefficients
    J_at_equilibrium = J.subs(equilibrium_point)

    a11 = J_at_equilibrium[0, 0]
    a12 = J_at_equilibrium[0, 1]
    a21 = J_at_equilibrium[1, 0]
    a22 = J_at_equilibrium[1, 1]

    # 6. The 'b' coefficients are zero for linearization at an equilibrium point
    b11 = 0
    b22 = 0
    
    # 7. Print the results as requested
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

if __name__ == '__main__':
    solve_linearization()