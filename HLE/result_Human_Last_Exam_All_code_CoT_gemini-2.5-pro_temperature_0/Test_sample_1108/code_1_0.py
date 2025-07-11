import sympy as sp

def solve_linearization():
    """
    This function calculates the linearization of the Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # Define symbolic variables
    S, F = sp.symbols('S F', real=True, positive=True)
    
    # Define the system of differential equations with given parameters
    # dS/dt = S(h - m*S/F)
    # dF/dt = F(a - b*F - c*S)
    # with a=2, b=c=h=m=1
    
    f = S * (1 - S / F)
    g = F * (2 - F - S)
    
    # Find the non-trivial equilibrium point (S_e, F_e) by solving f=0 and g=0
    # Since S > 0 and F > 0, the equations simplify to:
    # 1 - S/F = 0  => S = F
    # 2 - F - S = 0
    # Substituting S = F into the second equation gives 2 - 2F = 0, so F = 1.
    # Therefore, S = 1.
    # The non-trivial equilibrium point is (S_e, F_e) = (1, 1).
    
    Se, Fe = 1, 1
    
    # Compute the Jacobian matrix J = [[df/dS, df/dF], [dg/dS, dg/dF]]
    J = sp.Matrix([f, g]).jacobian(sp.Matrix([S, F]))
    
    # Evaluate the Jacobian at the equilibrium point (Se, Fe)
    J_eq = J.subs({S: Se, F: Fe})
    
    # The coefficients a_ij are the entries of the evaluated Jacobian matrix
    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]
    
    # For linearization at an equilibrium point, the constant terms b_ij are 0.
    b11 = 0
    b22 = 0
    
    # Print the results
    print(f"The non-trivial equilibrium point is (S, F) = ({Se}, {Fe}).")
    print("The coefficients of the linearized system are:")
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

solve_linearization()