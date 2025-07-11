import sympy

def solve_circle_problem():
    """
    Solves for r^2 based on the geometric constraints of two tangent circles.
    """
    # Define r as a positive symbolic variable
    r = sympy.Symbol('r', positive=True)
    # The radius of the second circle, D
    R = 2

    # Step 1: The two tangent lines are y = x + 1 and y = -x + 5.
    # They are perpendicular and intersect at the point P(2, 3).
    # The angle bisectors of these lines are x = 2 and y = 3.
    # The center of any circle tangent to both lines must lie on one of the bisectors.
    
    # Step 2: For a circle tangent to these specific lines, the distance (d) from its center
    # to the intersection point P(2,3) is related to its radius (rad) by the formula:
    # d = rad * sqrt(2).
    # So, for circle C: dist(C, P) = r * sqrt(2)
    # And for circle D: dist(D, P) = R * sqrt(2) = 2 * sqrt(2)
    
    # Step 3: To ensure a unique solution as implied by the problem's phrasing,
    # we consider the most symmetric configuration: the centers of the two circles lie on different angle bisectors.
    # Let center C lie on the line x=2 and center D lie on the line y=3.
    # The distance between C and D can be found using the Pythagorean theorem,
    # with the legs being the distances of C and D from the intersection point P.
    
    # dist(C, D)^2 = (dist(C, P))^2 + (dist(D, P))^2
    # Substituting the expressions from Step 2:
    dist_CP_sq = (r * sympy.sqrt(2))**2
    dist_DP_sq = (R * sympy.sqrt(2))**2
    dist_CD_sq_lhs = dist_CP_sq + dist_DP_sq
    
    # Step 4: The two circles are tangent to each other. For external tangency, the distance
    # between their centers is the sum of their radii.
    # dist(C, D) = r + R = r + 2
    # Squaring this gives:
    dist_CD_sq_rhs = (r + R)**2

    # Step 5: Form an equation by setting the two expressions for dist(C, D)^2 equal.
    # 2*r^2 + 8 = (r + 2)^2
    equation = sympy.Eq(dist_CD_sq_lhs, dist_CD_sq_rhs)
    
    # Solve the equation for r
    # 2*r**2 + 8 = r**2 + 4*r + 4
    # r**2 - 4*r + 4 = 0
    # (r - 2)**2 = 0
    solutions = sympy.solve(equation, r)
    
    # Since r must be positive, we take the positive solution.
    final_r = solutions[0]
    
    # Step 6: Calculate the final answer, r^2.
    r_squared = final_r**2

    # Print the steps and the final answer
    print("Let r be the radius of the first circle and R=2 be the radius of the second circle.")
    print("The two tangent lines y=x+1 and y=-x+5 are perpendicular and meet at P(2,3).")
    print("The centers of circles tangent to these lines must lie on the angle bisectors x=2 or y=3.")
    print("\nThe distance from a center to P is given by d = radius * sqrt(2).")
    print("So, dist(C,P)^2 = (r*sqrt(2))^2 = 2*r^2")
    print("And, dist(D,P)^2 = (2*sqrt(2))^2 = 8")
    print("\nAssuming centers C and D are on different bisectors, the squared distance between them is:")
    print("dist(C,D)^2 = dist(C,P)^2 + dist(D,P)^2")
    print(f"Equation Part 1: dist(C,D)^2 = {2}*r**{2} + {8}")
    print("\nSince the circles are tangent, the distance between centers is also the sum of the radii:")
    print("dist(C,D) = r + R = r + 2")
    print(f"Equation Part 2: dist(C,D)^2 = (r + {2})**{2} = r**{2} + {4}*r + {4}")
    print("\nSetting the two expressions equal:")
    print(f"{2}*r**{2} + {8} = r**{2} + {4}*r + {4}")
    print("Solving for r gives r = 2.")
    print("\nThe question asks for the value of r^2.")
    print(f"r^2 = {final_r}**{2} = {r_squared}")

if __name__ == '__main__':
    solve_circle_problem()