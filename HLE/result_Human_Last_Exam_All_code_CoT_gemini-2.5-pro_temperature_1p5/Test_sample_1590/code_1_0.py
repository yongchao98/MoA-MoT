import sympy

def solve_for_sliding_angle():
    """
    This script symbolically derives the expression for the angle at which the rod begins to slide.
    """
    # Define the symbolic variables based on the problem description.
    # M = mass, g = gravity, mu = coefficient of friction, theta = angle.
    M, g, mu, theta = sympy.symbols('M g mu theta')

    print("Step 1: Define the components of the gravitational force (Mg) relative to the tilted rod.")
    # The component of gravity perpendicular to the rod. This is balanced by the Normal Force (N).
    mg_perp = M * g * sympy.cos(theta)
    # The component of gravity parallel to the rod. This is the force pulling the rod down its slope.
    mg_parallel = M * g * sympy.sin(theta)
    print(f"Force component perpendicular to rod (N): {mg_perp}")
    print(f"Force component parallel to rod (f_required): {mg_parallel}\n")

    print("Step 2: Apply the condition for sliding.")
    print("Sliding begins when the parallel force component equals the maximum static friction force (mu * N).\n")

    # The equation for the onset of sliding: f_required = mu * N
    # We substitute the expressions from Step 1.
    final_equation = sympy.Eq(mg_parallel, mu * mg_perp)

    print("Step 3: Formulate and simplify the final equation.")
    print("The final equation relating the forces is:")
    print(f"{final_equation.lhs} = {final_equation.rhs}")
    
    # To simplify, we can divide both sides by M*g, which cancels out.
    # sin(theta) = mu * cos(theta)
    # Dividing by cos(theta) gives tan(theta) = mu.
    simplified_equation = sympy.Eq(sympy.tan(theta), mu)
    print("\nThis simplifies to:")
    print(str(simplified_equation))
    
    print("\nNote: The lengths L and l, and the mass M, do not affect the angle at which sliding begins.")
    print("There are no numerical constants from the problem description in the final equation.\n")
    
    # Step 4: Solve for theta.
    solution = sympy.solve(simplified_equation, theta)
    
    print("Step 4: The final expression for the angle theta is:")
    print(f"theta = {solution[0]}")

# Execute the function to print the derivation.
solve_for_sliding_angle()