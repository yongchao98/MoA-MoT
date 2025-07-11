import sympy

def find_separatrix():
    """
    This function verifies the analytical solution for the separatrix of the given system of ODEs.
    """
    # Define symbolic variables
    d_var, u_var = sympy.symbols('d u')

    # Define the functions F and G from the system d'=F(d,u), u'=G(d,u)
    F = 2*d_var**2 + (-3*u_var + 5*u_var**2)*d_var - u_var**3*(1-u_var)
    G = (u_var-1)*u_var**2

    print("The system of differential equations is:")
    print("d'(t) = 2*d^2 + (-3*u + 5*u^2)*d - u^3*(1 - u)")
    print("u'(t) = (u - 1)*u^2")
    print("-" * 30)

    print("We are looking for a separatrix, which is a special trajectory of the system.")
    print("From analytical calculations, the saddle point is at (d, u) = (-1, 1).")
    print("The separatrix is an invariant curve passing through the saddle point.")
    print("Our candidate solution is d = -u^2.")
    print("-" * 30)
    
    # Proposed separatrix equation
    separatrix_candidate = -u_var**2
    
    print(f"To verify the candidate d = {sympy.pretty(separatrix_candidate)}, we check if it satisfies the trajectory equation dd/du = F/G.")

    # 1. Calculate the LHS: d/du of the candidate solution
    lhs = sympy.diff(separatrix_candidate, u_var)
    print(f"\nLHS = d/du({sympy.pretty(separatrix_candidate)}) = {sympy.pretty(lhs)}")

    # 2. Calculate the RHS: substitute d = -u^2 into F/G
    rhs = F.subs(d_var, separatrix_candidate) / G
    # Simplify the RHS
    rhs_simplified = sympy.simplify(rhs)
    
    print(f"\nRHS = F/G with d = {sympy.pretty(separatrix_candidate)}.")
    # print(f"RHS before simplification = ({sympy.pretty(F.subs(d_var, separatrix_candidate))}) / ({sympy.pretty(G)})")
    print(f"RHS after simplification = {sympy.pretty(rhs_simplified)}")

    # 3. Check if LHS and RHS are equal
    if sympy.simplify(lhs - rhs_simplified) == 0:
        print("\nConclusion: LHS is equal to RHS.")
        print("The equation d = -u^2 is a valid invariant manifold (trajectory).")
        print("Since it passes through the saddle point (-1, 1), it is a separatrix.")
        
        # Output the final equation and its numbers
        print("-" * 30)
        print("The final equation for the separatrix is d = -u^2.")
        print("This can be written as: d + u^2 = 0.")
        print("\nDisplaying the numbers in the final equation (form: a*d + b*u^c = k):")
        d_coeff = 1
        u_coeff = 1
        u_power = 2
        constant = 0
        print(f"The equation can be expressed as: {d_coeff}*d + {u_coeff}*u^{u_power} = {constant}")

    else:
        print("\nVerification failed. The candidate is not a solution.")

if __name__ == '__main__':
    find_separatrix()