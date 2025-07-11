import sympy as sp

def find_and_verify_separatrix():
    """
    This function verifies the equation of the separatrix for the given system
    of differential equations and prints the resulting equation.

    The system of differential equations is:
    d'(t) = 2*d(t)**2 + (-3*u(t) + 5*u(t)**2)*d(t) - u(t)**3*(1 - u(t))
    u'(t) = (u(t) - 1)*u(t)**2

    Analysis of the system reveals a saddle point at (d, u) = (-1, 1).
    The separatrix is an invariant manifold passing through this saddle point.
    Through analytical methods, the equation for the separatrix is hypothesized to be d = -u**2.
    This script will programmatically verify this hypothesis.
    """
    print("Verifying the hypothesized separatrix equation...")
    
    # Define symbolic variables
    u = sp.Symbol('u')
    
    # The hypothesized separatrix equation is d = -u^2.
    # We treat d as a function of u, which is in turn a function of t.
    d = -u**2

    # From the second DE, we have the expression for u'(t) in terms of u.
    u_prime_t = (u - 1) * u**2

    # The left-hand side (LHS) of the first DE is d'(t).
    # We can calculate it using the chain rule: d'(t) = (dd/du) * u'(t).
    dd_du = sp.diff(d, u)
    lhs_d_prime = dd_du * u_prime_t

    # The right-hand side (RHS) of the first DE is given in the problem statement.
    # We substitute our hypothesized expression for d into it.
    rhs_d_prime = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1 - u)

    # To verify our hypothesis, we check if LHS equals RHS.
    # This is equivalent to checking if (LHS - RHS) simplifies to zero.
    verification_expression = lhs_d_prime - rhs_d_prime
    simplified_verification = sp.simplify(verification_expression)

    if simplified_verification == 0:
        print("Verification successful. The expression (LHS - RHS) simplifies to 0.")
        
        # The separatrix equation is of the form d = a*u^2 + b*u + c.
        # For d = -u^2, the coefficients are a=-1, b=0, c=0.
        a = -1
        b = 0
        c = 0
        
        print("The equation for the separatrix is d(u) = -u^2.")
        print("Displaying the equation with its numerical coefficients:")
        print(f"d = ({a})*u**2 + ({b})*u + ({c})")
    else:
        print("Verification failed. The hypothesized separatrix is not correct.")
        print(f"The simplified (LHS - RHS) expression is: {simplified_verification}")

# Execute the function
find_and_verify_separatrix()
<<<d = -u**2>>>