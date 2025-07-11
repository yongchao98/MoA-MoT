import sympy

def solve_particle_emitter_problem():
    """
    Calculates the minimum ratio of the cube of the surface area to the square of the volume
    for the region traversed by particles from an emitter.
    """
    # Step 1 & 2: Define symbols and geometric properties
    # x is the dimensionless ratio v^2 / (2gh), where v is initial speed,
    # g is gravity, and h is emitter height.
    x = sympy.Symbol('x', positive=True)
    pi = sympy.pi

    # We can set the initial height h=1 without loss of generality, as it will
    # cancel out from the final dimensionless ratio A^3 / V^2.
    h = 1

    # H is the total height of the "parabola of safety" from its vertex to the ground plane.
    # k is the curvature parameter of the paraboloid (its equation is z = H - k*r^2).
    # R2 is the square of the radius of the base where the paraboloid intersects the ground.
    H = h * (1 + x)
    k = 1 / (4 * h * x)
    R2 = H / k

    # Step 3: Define the Volume (V) and Surface Area (A)
    # The volume of a paraboloid segment is (1/2) * base_area * height.
    V = (pi / 2) * R2 * H

    # The total surface area is the sum of the circular base and the curved paraboloid surface.
    A_base = pi * R2
    # The standard formula for the surface area of a paraboloid of revolution is used for the curved part.
    A_paraboloid = (pi / (6 * k**2)) * ((1 + 4 * k**2 * R2)**(sympy.S(3)/2) - 1)
    A = A_base + A_paraboloid

    # Step 4: Form the ratio and simplify it
    # We are asked to find the minimum of A^3 / V^2.
    Ratio = (A**3) / (V**2)
    # Simplifying this expression shows it is independent of the initial height 'h'.
    Ratio = sympy.simplify(Ratio)

    # Step 5: Find the value of x that minimizes the ratio
    # To find the minimum, we take the derivative with respect to x and set it to zero.
    diff_Ratio = sympy.diff(Ratio, x)
    
    # The derivative is a complex fraction. The minimum occurs when its numerator is zero.
    num, den = sympy.fraction(sympy.simplify(diff_Ratio))
    
    # We solve the equation numerator = 0 for x.
    solutions = sympy.solve(num, x)
    
    # We select the single physically meaningful solution (positive and real).
    x_min = [sol for sol in solutions if sol.is_real and sol.is_positive][0]

    # Step 6: Substitute the optimal x back into the ratio to find the minimum value
    min_Ratio = Ratio.subs(x, x_min)
    final_answer = sympy.simplify(min_Ratio)

    # Step 7: Display the final result and its components as requested
    print("The minimum ratio of the cube of the surface area to the square of the volume has been determined.")
    print(f"The exact symbolic expression for the minimum ratio is: {final_answer}")
    
    # The final answer is of the form (A + B*sqrt(C))*pi.
    # We expand the expression to find the integer components A, B, and C.
    expanded_answer = sympy.expand(final_answer / pi)
    
    # Extract coefficients programmatically
    # The expression is of the form const + term_with_sqrt
    const_part = expanded_answer.as_ordered_terms()[0]
    sqrt_part = expanded_answer.as_ordered_terms()[1]
    
    A_val = const_part
    B_val = sqrt_part.as_coeff_mul()[0]
    C_val = sqrt_part.as_coeff_mul()[1].args[0]

    print("\nThe final answer is in the form (A + B*sqrt(C))*pi.")
    print("As requested, the integer components of this equation are:")
    print(f"A = {A_val}")
    print(f"B = {B_val}")
    print(f"C = {C_val}")
    
    numerical_answer = final_answer.evalf()
    print(f"\nThe numerical value of this minimum ratio is approximately: {numerical_answer:.5f}")

if __name__ == '__main__':
    solve_particle_emitter_problem()