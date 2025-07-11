import sympy

def solve_linearization():
    """
    This function calculates the linearization of the given predator-prey system
    at its non-trivial equilibrium point.
    """
    # Step 1: Define symbolic variables and parameters
    S, F = sympy.symbols('S F')
    a, b, c, h, m = (2, 1, 1, 1, 1)

    # Define the system of differential equations
    dS_dt = S * (h - m * S / F)
    dF_dt = F * (a - b * F - c * S)

    # Step 2: Find the non-trivial equilibrium point
    # We solve dS/dt = 0 and dF/dt = 0
    # The non-trivial solution requires S > 0 and F > 0.
    # From dS/dt = 0, we get h - m*S/F = 0 => S = (h/m)*F. With our params, S = F.
    # Substitute S = F into dF/dt = 0: F*(a - b*F - c*F) = 0
    # For F > 0, we have a - (b+c)*F = 0 => F = a/(b+c).
    # With our params, F = 2/(1+1) = 1.
    # So, the equilibrium point is S=1, F=1.
    S_eq, F_eq = 1, 1

    # Step 3: Compute the Jacobian Matrix of the system
    # The functions are f(S,F) = dS/dt and g(S,F) = dF/dt
    system_matrix = sympy.Matrix([dS_dt, dF_dt])
    jacobian_matrix = system_matrix.jacobian([S, F])
    
    # Step 4: Evaluate the Jacobian matrix at the equilibrium point
    J_eq = jacobian_matrix.subs({S: S_eq, F: F_eq})

    a11 = J_eq[0, 0]
    a12 = J_eq[0, 1]
    a21 = J_eq[1, 0]
    a22 = J_eq[1, 1]
    
    # Step 5: Determine the constant terms b_ij
    # Since we linearize at an equilibrium point, the system's value is (0,0).
    b11 = dS_dt.subs({S: S_eq, F: F_eq})
    b22 = dF_dt.subs({S: S_eq, F: F_eq})

    # Step 6: Print the results
    print("The coefficients for the linearized system are:")
    print(f"a11 = {a11}")
    print(f"a12 = {a12}")
    print(f"a21 = {a21}")
    print(f"a22 = {a22}")
    print(f"b11 = {b11}")
    print(f"b22 = {b22}")
    
    print("\nThe linearized system in matrix form is:")
    print(f"| x'(t) | = | {a11}  {a12} | | x(t) | + | {b11} |")
    print(f"| y'(t) |   | {a21} {a22} | | y(t) |   | {b22} |")
    
    print("\nThe linearized equations are:")
    print(f"x'(t) = ({a11})*x(t) + ({a12})*y(t) + {b11}")
    print(f"y'(t) = ({a21})*x(t) + ({a22})*y(t) + {b22}")
    
    # Return the answer string for the final block
    return f"a11 = {a11}, a12 = {a12}, a21 = {a21}, a22 = {a22}, b11 = {b11}, b22 = {b22}"

# Execute the function and capture the result for the final answer
final_answer = solve_linearization()
print(f"\n<<<{final_answer}>>>")
