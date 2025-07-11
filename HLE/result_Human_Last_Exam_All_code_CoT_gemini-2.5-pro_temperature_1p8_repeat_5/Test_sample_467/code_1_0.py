import sympy

def solve_minimal_surface_index():
    """
    Calculates the Morse index of a minimal surface M with Gauss map g(z) = z/(z^3+2).

    The Morse index is calculated using the Nayatani formula:
    Index = I(g) + k - 1
    where:
    k = number of ends of the surface
    I(g) = total number of poles of the 1-form dg/g on the Riemann sphere
    """
    z = sympy.Symbol('z')

    # Define the Gauss map g(z)
    g = z / (z**3 + 2)
    print(f"The Gauss map is g(z) = {g}\n")

    # Step 1: Calculate k, the number of ends.
    # Ends are the poles of the height differential, which occur where g=0 or g=infinity
    # on the Riemann sphere (C U {infinity}).
    
    # Poles of g are the roots of the denominator.
    g_num, g_den = g.as_numer_denom()
    poles_of_g_in_C = sympy.solve(g_den, z)
    num_poles_g_in_C = len(poles_of_g_in_C)
    print(f"The denominator of g(z) is: {g_den}")
    print(f"The roots of the denominator are the solutions to z^3 + 2 = 0.")
    print(f"Number of finite poles of g(z): {num_poles_g_in_C}")

    # Poles of 1/g are the roots of the numerator.
    poles_of_inv_g_in_C = sympy.solve(g_num, z)
    num_poles_inv_g_in_C = len(poles_of_inv_g_in_C)
    print(f"\nThe numerator of g(z) is: {g_num}")
    print(f"The root is z = 0.")
    print(f"Number of finite poles of 1/g(z): {num_poles_inv_g_in_C}")
    
    # Analyze behavior at infinity by substituting z = 1/w and checking w=0.
    w = sympy.Symbol('w')
    g_at_inf = g.subs(z, 1/w)
    
    is_pole_g_at_inf = sympy.limit(abs(g_at_inf), w, 0) == sympy.oo
    is_pole_inv_g_at_inf = sympy.limit(abs(1/g_at_inf), w, 0) == sympy.oo

    # The point at infinity is an end if g=0 or g=oo there.
    # Since limit(g_at_inf, w, 0) is 0, it's a pole of 1/g, hence an end.
    num_ends_at_inf = 1 if is_pole_g_at_inf or is_pole_inv_g_at_inf else 0
    print(f"\nAnalyzing the point at infinity:")
    print(f"g(1/w) = {sympy.simplify(g_at_inf)}")
    print(f"As z -> oo (w -> 0), g(z) -> {sympy.limit(g_at_inf, w, 0)}.")
    print("Since g(z) -> 0, infinity is a pole of 1/g, and thus is an end.")
    print(f"Number of ends at infinity: {num_ends_at_inf}")

    k = num_poles_g_in_C + num_poles_inv_g_in_C + num_ends_at_inf
    print(f"\nTotal number of ends k = {num_poles_g_in_C} + {num_poles_inv_g_in_C} + {num_ends_at_inf} = {k}\n")
    
    # Step 2: Calculate I(g), the number of poles of dg/g.
    # This is the total number of poles of the function g'(z)/g(z).
    dg_dz = sympy.diff(g, z)
    omega_func = sympy.simplify(dg_dz / g)
    print(f"The logarithmic derivative is g'(z)/g(z) = {omega_func}")

    # Find finite poles of omega_func
    omega_num, omega_den = omega_func.as_numer_denom()
    num_poles_omega_in_C = sympy.degree(omega_den, z)
    print(f"\nThe denominator of g'(z)/g(z) is: {omega_den}")
    print(f"The degree of the denominator is {num_poles_omega_in_C}, which gives the number of finite poles.")

    # Check for pole at infinity
    g_inf = g.subs(z, 1/w)
    dg_inf_dw = sympy.diff(g_inf, w)
    omega_func_inf = sympy.simplify(dg_inf_dw / g_inf)
    is_pole_omega_at_inf = sympy.limit(abs(omega_func_inf), w, 0) == sympy.oo
    num_poles_omega_at_inf = 1 if is_pole_omega_at_inf else 0
    print(f"\nAnalyzing g'(z)/g(z) at infinity:")
    # We transform the function f(z) = g'(z)/g(z) to a function of w
    # Note that g'(z) = g'(w) * dw/dz = g'(w) * (-w^2). This factor does not affect pole existence.
    print(f"The corresponding function in w is: g'(w)/g(w) = {omega_func_inf}")
    print(f"As w -> 0, this function tends to oo, so there is a pole at infinity.")
    print(f"Number of poles of dg/g at infinity: {num_poles_omega_at_inf}")

    I_g = num_poles_omega_in_C + num_poles_omega_at_inf
    print(f"\nTotal number of poles of dg/g, I(g) = {num_poles_omega_in_C} + {num_poles_omega_at_inf} = {I_g}\n")

    # Step 3: Calculate the Morse Index
    index = I_g + k - 1
    print("Final Calculation using Nayatani's formula:")
    print(f"Index = I(g) + k - 1")
    print(f"Index = {I_g} + {k} - 1 = {index}")
    
solve_minimal_surface_index()