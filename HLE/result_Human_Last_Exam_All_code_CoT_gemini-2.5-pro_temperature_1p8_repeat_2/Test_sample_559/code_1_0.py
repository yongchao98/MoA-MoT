import sympy as sp

def find_separatrix():
    """
    This function finds the separatrix for the given system of ODEs by hypothesizing
    a polynomial relationship between d and u and solving for its coefficients.
    """
    # Define u and c as symbolic variables
    u, c = sp.symbols('u c')

    # The given system of differential equations is:
    # d'(t) = 2*d(t)^2 + (-3*u(t) + 5*u(t)^2)*d(t) - u(t)^3*(1-u(t))
    # u'(t) = u(t)^2 * (u(t) - 1)

    print("Step 1: Formulate a hypothesis for the separatrix.")
    print("We are looking for a special solution (an invariant manifold) to the system.")
    print("Let's hypothesize a solution of the form: d = c * u**2.")
    d_hypothesis = c * u**2

    print("\nStep 2: Substitute the hypothesis into the system of equations.")
    print("First, express d'(t) in terms of u(t) and u'(t) using the chain rule:")
    # d/dt = (d(d)/du) * (du/dt) = (2*c*u) * u_prime
    # We substitute the given expression for u'(t) = u^2 * (u - 1)
    u_prime_expr = u**2 * (u - 1)
    d_prime_from_hypothesis = (2 * c * u) * u_prime_expr
    print(f"d'(t) = 2*c*u * u'(t) = 2*c*u * ({sp.simplify(u_prime_expr)}) = {sp.simplify(d_prime_from_hypothesis)}")

    print("\nNext, substitute d = c * u**2 into the right-hand side (RHS) of the first ODE:")
    rhs_ode = 2 * d_hypothesis**2 + (-3 * u + 5 * u**2) * d_hypothesis - u**3 * (1 - u)
    print(f"RHS = {sp.simplify(rhs_ode)}")

    print("\nStep 3: Equate the two expressions and solve for the constant c.")
    # For d = c*u**2 to be a solution, d'(t) must equal the RHS for all u.
    # So, d_prime_from_hypothesis - rhs_ode = 0.
    equation = sp.simplify(d_prime_from_hypothesis - rhs_ode)
    print("The equation (d'(t) - RHS) = 0 simplifies to:")
    print(f"{equation} = 0")
    
    # This polynomial in u must be zero for all u. This means each of its coefficients must be zero.
    # The equation is of the form: A*u**4 + B*u**3 = 0.
    # So, A=0 and B=0.
    # A = -2*c**2 - 3*c - 1
    # B = c + 1
    # From B=0, we get c = -1. We verify this in A.
    # -2*(-1)**2 - 3*(-1) - 1 = -2 + 3 - 1 = 0.
    # So, c=-1 is the correct value.

    # We can also solve programmatically.
    solution = sp.solve(equation, c)
    if not solution:
        print("No constant c found. The hypothesis is incorrect.")
        return
        
    c_val = solution[0]
    print(f"\nSolving for c, we find c = {c_val}.")

    print("\nStep 4: Write the final equation for the separatrix.")
    final_equation = d_hypothesis.subs(c, c_val)
    coeff, u_base, u_exp = final_equation.as_three_terms()
    
    print("Substituting c back into the hypothesis d = c * u**2 gives the separatrix equation.")
    print(f"Final Equation: d = {int(coeff)} * u**{u_exp}")

    print("\nThe numbers in the final equation are:")
    print(f"Coefficient: {int(coeff)}")
    print(f"Exponent: {u_exp}")


if __name__ == '__main__':
    find_separatrix()
    final_equation_str = "d = -u**2"
    print(f"\n<<<{final_equation_str}>>>")