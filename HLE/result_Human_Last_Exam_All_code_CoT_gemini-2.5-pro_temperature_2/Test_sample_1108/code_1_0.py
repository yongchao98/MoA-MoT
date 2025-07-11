import sympy

def solve_leslie_linearization():
    """
    Calculates the coefficients of the linearized Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # 1. Define the system with given parameters
    S, F = sympy.symbols('S F', real=True, positive=True)
    a, b, c, h, m = 2, 1, 1, 1, 1
    
    # Define the functions for dS/dt and dF/dt
    f_S = S * (h - m * S / F)
    f_F = F * (a - b * F - c * S)

    # 2. Find the non-trivial equilibrium point
    # We solve f_S = 0 and f_F = 0 for S > 0, F > 0
    # From f_S = 0 => S/F = 1 => S = F
    # Substitute S = F into f_F = 0 => F(a - bF - cF) = 0 => a = (b+c)F
    # F_e = a / (b+c) = 2 / (1+1) = 1
    # S_e = F_e = 1
    S_e, F_e = 1, 1

    # 3. Compute the Jacobian Matrix
    system_matrix = sympy.Matrix([f_S, f_F])
    variables_matrix = sympy.Matrix([S, F])
    jacobian_matrix = system_matrix.jacobian(variables_matrix)

    # 4. Evaluate the Jacobian at the equilibrium point (S_e, F_e)
    jacobian_at_equilibrium = jacobian_matrix.subs({S: S_e, F: F_e})
    
    a11 = jacobian_at_equilibrium[0, 0]
    a12 = jacobian_at_equilibrium[0, 1]
    a21 = jacobian_at_equilibrium[1, 0]
    a22 = jacobian_at_equilibrium[1, 1]

    # 5. Determine the constant terms b11 and b22
    # These are the values of the functions at the equilibrium, which are zero.
    b11 = f_S.subs({S: S_e, F: F_e})
    b22 = f_F.subs({S: S_e, F: F_e})

    # Print the coefficients of the final equation
    print("The non-trivial equilibrium point (S_e, F_e) is ({}, {}).".format(S_e, F_e))
    print("\nThe linearized system is given by:")
    print("x'(t) = a11*x(t) + a12*y(t) + b11")
    print("y'(t) = a21*x(t) + a22*y(t) + b22")
    print("\nThe coefficients are:")
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

if __name__ == '__main__':
    solve_leslie_linearization()
    # For the final answer format
    print("\n<<<-1, 1, -1, -1, 0, 0>>>")
