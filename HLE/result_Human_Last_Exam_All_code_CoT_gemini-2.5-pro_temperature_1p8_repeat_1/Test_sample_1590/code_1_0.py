import sympy as sp

def solve_sliding_angle():
    """
    This function symbolically derives the expression for the angle at which a
    tilted rod begins to slide off a frictional corner.
    """
    # 1. Define the symbolic variables
    # M: mass of the rod
    # g: acceleration due to gravity
    # N: normal force from the corner on the rod
    # f: friction force from the corner on the rod
    # mu: coefficient of static friction
    # theta: angle of the rod with the horizontal
    M, g, N, f, mu, theta = sp.symbols('M g N f mu theta', real=True, positive=True)

    print("Derivation for the angle at which the rod begins to slide:")
    print("----------------------------------------------------------\n")

    # 2. Force balance equations based on a tilted coordinate system (parallel and perpendicular to the rod)
    print("1. Resolve the gravitational force Mg into components parallel and perpendicular to the rod:")
    # The component of gravity perpendicular to the rod is balanced by the Normal force N.
    eq_normal_force = sp.Eq(N, M * g * sp.cos(theta))
    # The component of gravity parallel to the rod is balanced by the Friction force f.
    eq_friction_force = sp.Eq(f, M * g * sp.sin(theta))
    
    print(f"   - Perpendicular component (balanced by Normal Force N): N = M * g * cos(theta)")
    print(f"   - Parallel component (balanced by Friction Force f): f = M * g * sin(theta)\n")
    
    # 3. State the condition for the onset of sliding
    # Sliding begins when the static friction force f reaches its maximum value, f_max = mu * N.
    eq_sliding_condition = sp.Eq(f, mu * N)
    print("2. The rod begins to slide when the static friction force equals its maximum value:")
    print(f"   - Sliding condition: {eq_sliding_condition}\n")

    # 4. Substitute and solve for theta
    # Substitute the expressions for f and N into the sliding condition.
    final_equation = eq_sliding_condition.subs([
        (N, eq_normal_force.rhs),
        (f, eq_friction_force.rhs)
    ])
    
    print("3. Substitute the expressions for N and f into the sliding condition:")
    print(f"   Equation: {final_equation.lhs} = {final_equation.rhs}\n")
    
    print("4. The terms for mass 'M' and gravity 'g' cancel out, simplifying the equation to:")
    print(f"   sin(theta) = mu * cos(theta)\n")
    
    print("5. To find the final expression for the angle, we rearrange the equation:")
    print("   sin(theta) / cos(theta) = mu")
    print("   tan(theta) = mu\n")
    
    # The final expression for theta
    final_expression = sp.atan(mu)
    
    print("Final equation for the angle theta:")
    # We print each part of the final equation as requested.
    part1 = "theta"
    part2 = " = "
    part3 = "arctan"
    part4 = "("
    part5 = "mu"
    part6 = ")"
    print(f"{part1}{part2}{part3}{part4}{part5}{part6}")

if __name__ == "__main__":
    solve_sliding_angle()