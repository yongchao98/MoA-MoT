import sympy

def solve_for_force():
    """
    This script symbolically derives the force per unit area on the conducting plane.
    It follows the steps outlined in the plan:
    1. Define potentials in each region.
    2. Apply boundary conditions to solve for unknown coefficients.
    3. Calculate the magnetic field at the conductor surface.
    4. Calculate the force from the magnetic pressure.
    """

    # 1. Define symbols and general solutions for the magnetic vector potential A_z
    x, y, a, d = sympy.symbols('x y a d', real=True, positive=True)
    K0, mu0, mu = sympy.symbols('K_0, mu_0, mu', real=True, positive=True)
    
    # Unknown coefficients for the potentials
    C, D1 = sympy.symbols('C D1')
    
    # Potential in Region 1 (0 < x < d, air gap).
    # The form sinh(a*(x-d)) is chosen to automatically satisfy the boundary condition
    # B_x = dA_z/dy = 0 at x=d.
    A_z1 = C * sympy.sinh(a * (x - d)) * sympy.sin(a * y)
    
    # Potential in Region 2 (x < 0, magnetic material).
    # The form exp(a*x) is chosen to ensure the field does not diverge at x -> -oo.
    A_z2 = D1 * sympy.exp(a * x) * sympy.sin(a * y)
    
    # 2. Apply boundary conditions at the x=0 interface
    
    # B = curl(A) -> Bx = diff(Az, y), By = -diff(Az, x)
    B_x1 = sympy.diff(A_z1, y)
    B_y1 = -sympy.diff(A_z1, x)
    H_y1 = B_y1 / mu0
    
    B_x2 = sympy.diff(A_z2, y)
    B_y2 = -sympy.diff(A_z2, x)
    H_y2 = B_y2 / mu
    
    # BC 1: Normal component of B is continuous: B_x1(0) = B_x2(0)
    eq1 = sympy.Eq(B_x1.subs(x, 0), B_x2.subs(x, 0))
    # This simplifies to C * sinh(-ad) = D1, or D1 = -C * sinh(ad)
    # Let's solve it formally
    sol_D1 = sympy.solve(eq1, D1)[0]

    # BC 2: Tangential component of H is discontinuous by K: H_y1(0) - H_y2(0) = K0 * sin(a*y)
    eq2 = sympy.Eq(H_y1.subs(x, 0) - H_y2.subs(x, 0), K0 * sympy.sin(a * y))
    
    # Substitute D1 into eq2 and solve for C
    eq2_sub = eq2.subs(D1, sol_D1)
    sol_C = sympy.solve(eq2_sub, C)[0]

    # 3. Calculate the magnetic field B at the conductor surface (x=d)
    # At x=d, B_x is zero. We only need B_y.
    B_y_at_d = B_y1.subs(x, d)

    # Substitute the solved coefficient C back into the expression for B_y
    B_y_final = B_y_at_d.subs(C, sol_C)
    
    # The total magnetic field magnitude squared at x=d is just B_y_final^2
    B_sq = sympy.simplify(B_y_final**2)
    
    # 4. Calculate the force per unit area using magnetic pressure P_m = B^2 / (2*mu_0)
    # The force is in the +x direction, pushing the conductor.
    force_per_area_magnitude = B_sq / (2 * mu0)
    
    final_force_expression = sympy.simplify(force_per_area_magnitude)
    
    # 5. Print the result in a format similar to the multiple-choice options.
    # We will build the string to represent the final formula clearly.
    numerator_str = f"μ₀ * K₀**2 * sin(a*y)**2"
    denominator_str = f"(cosh(a*d) + (μ₀/μ) * sinh(a*d))**2"
    full_expression = f"({numerator_str}) / (2 * {denominator_str}) * î_x"

    print("Derived Force Per Unit Area:")
    sympy.pprint(final_force_expression)
    print("\nFormatted Expression:")
    print(f"f/area = {full_expression}")

    # The result corresponds to option C.
    print("\nComparing the derived expression with the answer choices, the correct choice is C.")

solve_for_force()