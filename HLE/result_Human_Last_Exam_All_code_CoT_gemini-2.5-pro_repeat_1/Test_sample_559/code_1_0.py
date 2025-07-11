import sympy as sp

def solve_separatrix():
    """
    This function verifies the proposed separatrix for the given system of ODEs
    and prints the final equation.
    """
    # Define symbolic variables
    # We use u as a symbol, not a function of t, for simplicity in phase plane analysis.
    # d is treated as a function of u.
    u = sp.Symbol('u')
    d = sp.Function('d')(u)

    # Define the system of ODEs in phase plane form: dd/du = d_prime / u_prime
    # d_prime = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    # u_prime = u**2 * (u-1)
    # So, dd/du = (2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)) / (u**2 * (u-1))

    # We propose the separatrix is given by the equation d = -u**2.
    # Let's verify if this is a solution to the phase plane ODE.
    
    proposed_d = -u**2
    
    # Calculate the left-hand side (LHS): dd/du
    lhs = sp.diff(proposed_d, u)

    # Calculate the right-hand side (RHS) by substituting d = -u**2 into the expression for dd/du
    d_prime_subs = 2*proposed_d**2 + (-3*u + 5*u**2)*proposed_d - u**3*(1-u)
    u_prime = u**2 * (u - 1)
    
    # We need to handle the case u=1 and u=0 where u_prime is zero.
    # The simplification will work for u not equal to 0 or 1.
    rhs = sp.simplify(d_prime_subs / u_prime)

    print("Verifying the proposed separatrix: d = -u**2")
    print("-" * 40)
    print(f"The derivative dd/du from the proposed solution is: {lhs}")
    print(f"The expression for dd/du from the ODEs is: {rhs}")
    
    # Check if the LHS and RHS are equal
    if sp.simplify(lhs - rhs) == 0:
        print("\nThe verification is successful: LHS equals RHS.")
        print("This confirms that d = -u**2 is a trajectory of the system.")
        
        # An equilibrium point is a saddle if its linearization has eigenvalues with opposite signs.
        # This analysis (done manually) shows the point (-1, 1) is a saddle.
        # Now, we check if our separatrix passes through the saddle point (-1, 1).
        u_saddle, d_saddle = 1, -1
        if proposed_d.subs(u, u_saddle) == d_saddle:
            print(f"The curve d = -u**2 passes through the saddle point ({d_saddle}, {u_saddle}).")
            print("Therefore, it is a separatrix.")
            
            # Print the final equation with its coefficients as requested.
            # The equation is d = -u**2, which can be written as d + u**2 = 0.
            d_coeff = 1
            u_sq_coeff = 1
            exponent = 2
            constant = 0
            
            print("\nThe equation for the separatrix is:")
            print(f"d = -u**{exponent}")
            print("In the standard form Ax + By = C, it can be written as:")
            print(f"({d_coeff})*d + ({u_sq_coeff})*u**{exponent} = {constant}")

        else:
            print("The curve does not pass through the saddle point.")
    else:
        print("\nVerification failed. The proposed solution is not a separatrix.")

if __name__ == '__main__':
    solve_separatrix()