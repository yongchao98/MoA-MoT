import sympy as sp

def solve_linearization():
    """
    This function performs the linearization of the given predator-prey system.
    1. It defines the system equations and parameters.
    2. It calculates the non-trivial equilibrium point.
    3. It computes the Jacobian matrix (partial derivatives).
    4. It evaluates the Jacobian at the equilibrium point to find the coefficients a_ij.
    5. It identifies the coefficients b_ij, which are zero for linearization at an equilibrium.
    6. It prints the coefficients and the final linearized equations.
    """
    
    # Define parameters and symbols
    S, F = sp.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1
    
    # Define the system equations
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)
    
    # Solve for the non-trivial equilibrium point (S > 0, F > 0)
    # S(h - mS/F) = 0 => h - mS/F = 0 => hF = mS => S = (h/m)F = F
    # F(a - bF - cS) = 0 => a - bF - cS = 0
    # Substituting S = F into the second equation:
    # a - bF - cF = 0 => a = (b+c)F => F = a / (b+c)
    Fe = a / (b + c)
    Se = Fe
    
    # Compute the Jacobian matrix elements (a_ij)
    a11_expr = sp.diff(dS_dt, S)
    a12_expr = sp.diff(dS_dt, F)
    a21_expr = sp.diff(dF_dt, S)
    a22_expr = sp.diff(dF_dt, F)
    
    # Evaluate the Jacobian at the equilibrium point (Se, Fe)
    a11 = a11_expr.subs([(S, Se), (F, Fe)])
    a12 = a12_expr.subs([(S, Se), (F, Fe)])
    a21 = a21_expr.subs([(S, Se), (F, Fe)])
    a22 = a22_expr.subs([(S, Se), (F, Fe)])
    
    # The b_ij terms are the values of the functions at equilibrium, which are zero.
    b11 = 0
    b22 = 0
    
    # Print the coefficients
    print("The coefficients for the linearized system are:")
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")

    # Print the final equations with the numbers
    print("\nThe linearized equations are:")
    print(f"x'(t) = ({a11}) * x(t) + ({a12}) * y(t) + ({b11})")
    print(f"y'(t) = ({a21}) * x(t) + ({a22}) * y(t) + ({b22})")

# Run the function
solve_linearization()
print("<<<-1, 1, -1, -1, 0, 0>>>")