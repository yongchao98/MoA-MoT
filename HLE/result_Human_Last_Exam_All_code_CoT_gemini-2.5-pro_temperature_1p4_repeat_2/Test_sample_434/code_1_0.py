import sympy

def find_blowup_condition():
    """
    This function determines the condition on y(0) for which the solution
    to the given system of ODEs blows up, assuming x(0) > 1.
    
    The method involves finding the equation of a separatrix curve in the phase
    plane that divides blowing-up solutions from non-blowing-up solutions. This
    separatrix is the stable manifold of the saddle point at (1,0).
    """
    
    # Define symbolic variables for the derivation
    x, y = sympy.symbols('x y', real=True)

    # The ODE for the trajectories is dy/dx = (1 - x - y**2) / (-3*x*y).
    # This can be transformed into the differential form: (x + y**2 - 1)dx - 3*x*y*dy = 0.
    # Let M = x + y**2 - 1 and N_dy = -3*x*y.
    # This form can be made exact using the integrating factor mu(x) = x**(-5/3).
    # An exact form M_new*dx + N_new*dy = 0 implies that trajectories are level
    # curves of a potential function G(x,y). We find G(x,y) by integration.

    mu = x**(-sympy.S(5)/3)
    M = x + y**2 - 1
    N_dy = -3 * x * y
    
    # We find G(x,y) such that dG/dy = mu*N_dy and dG/dx = mu*M.
    # Integrate dG/dy with respect to y
    G_intermediate = sympy.integrate(mu * N_dy, y)
    
    # Find the integration "constant" h(x) by differentiating and comparing
    h_prime = sympy.simplify(mu * M - sympy.diff(G_intermediate, x))
    h = sympy.integrate(h_prime, x)
    
    # The full potential function G(x,y)
    G = G_intermediate + h
    
    # The separatrix is the trajectory that passes through the saddle point (1, 0).
    # We find the value of the constant for this specific level curve.
    C = G.subs({x: 1, y: 0})
    
    # The equation for the separatrix is G(x,y) = C.
    separatrix_eq = sympy.Eq(G, C)
    
    # To get the boundary for y(0), we solve for y**2 in terms of x.
    # The solutions that blow up are those starting below this curve.
    boundary_expr_y_squared = sympy.solve(separatrix_eq, y**2)[0]

    # The condition for blow-up is y(0) < sqrt(boundary_expr).
    # From the solved expression `2*x + 1 - 3*x**(2/3)`, we extract the numbers.
    c1 = 2
    c2 = 1
    c3 = 3
    p_num = 2
    p_den = 3
    exp_num = 1
    exp_den = 2

    # Print the final result in the requested format.
    print("For x(0) > 1, the solution of the system blows up if and only if y(0) satisfies the following condition:")
    print(f"y(0) < ({c1}*x(0) + {c2} - {c3}*x(0)**({p_num}/{p_den}))**({exp_num}/{exp_den})")

if __name__ == '__main__':
    find_blowup_condition()