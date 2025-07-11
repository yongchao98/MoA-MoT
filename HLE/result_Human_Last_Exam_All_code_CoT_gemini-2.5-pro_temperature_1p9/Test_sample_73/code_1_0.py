import sympy
import numpy as np

def solve_pentagon_genus():
    """
    Calculates the genus of the configuration space of a hinged equilateral pentagon
    with one side fixed in the plane.
    """
    print("Step-by-step derivation of the genus of the pentagon's configuration space:")
    print("-" * 70)

    # Step 1 & 2: Define the system and parameterization
    print("Step 1: The configuration is determined by the positions of three free vertices.")
    print("Step 2: We parameterize the system using the angles of two links.")
    print("Let the fixed side be V1=(0,0) and V2=(L,0).")
    print("The position of V3 depends on an angle theta_2. The position of V5 on theta_1.")
    print("The space of these parameters (theta_1, theta_2) is a torus T^2.")
    print("-" * 70)

    # Step 3: Define the constraint for the existence of V4
    print("Step 3: A valid position for V4 exists only if dist(V3, V5) <= 2L.")
    print("This condition defines a region S on the torus.")
    print("The boundary of S is where dist(V3, V5)^2 = (2L)^2, which simplifies to the equation:")
    print("2*cos(theta_2) - 2*cos(theta_1) - 2*cos(theta_1 - theta_2) = 1")
    print("-" * 70)

    # Step 4: The configuration space as a double
    print("Step 4: The full configuration space M is the 'double' of the region S.")
    print("Its genus g(M) is given by the formula: g(M) = 2*g(S) + k - 1")
    print("where g(S) is the genus of S and k is the number of its boundary components.")
    print("-" * 70)

    # Step 5: Analyze the boundary equation to find k
    print("Step 5: We analyze the boundary to find k.")
    print("The boundary equation can be solved for theta_2 for a given theta_1.")
    print("This has solutions only if a certain condition on theta_1 holds.")
    print("This condition reduces to the quadratic inequality: 4*cos^2(theta_1) + 12*cos(theta_1) - 7 <= 0.")
    
    x = sympy.Symbol('x')
    equation = 4*x**2 + 12*x - 7
    roots = sympy.solve(equation, x)
    
    print(f"\nLet x = cos(theta_1). The equation is {equation} = 0.")
    print(f"Using sympy to find the roots of the quadratic equation: {roots}")
    
    # Sympy gives exact solutions, let's get floats for interpretation.
    root1_float = float(roots[0])
    root2_float = float(roots[1])

    print(f"The roots are {root1_float:.2f} and {root2_float:.2f}. So the inequality holds for x between them.")
    print(f"Since -1 <= cos(theta_1) <= 1, the valid range is -1 <= cos(theta_1) <= 0.5.")
    print("This corresponds to theta_1 in the interval [pi/3, 5pi/3].")
    print("\nBecause theta_1 is restricted to an interval of width less than 2*pi, the boundary")
    print("curve is a single, non-self-intersecting loop. Therefore, k = 1.")
    print("-" * 70)

    # Step 6: Determine the topology of S to find g(S)
    print("Step 6: We determine the genus g(S).")
    print("A single contractible loop separates the torus into two regions:")
    print("  1. A region homeomorphic to a disk (genus 0).")
    print("  2. A region homeomorphic to a punctured torus (genus 1).")
    print("\nWe check which region S is by testing a point. Let's test (theta_1, theta_2) = (0,0).")
    print("The inequality is 2*cos(0) - 2*cos(0) - 2*cos(0) - 1 <= 0, which is -1 <= 0. True.")
    print("The point (0,0) is in S. The boundary is restricted to theta_1 in [pi/3, 5pi/3].")
    print("Since the point (0,0) is outside the region enclosed by the boundary curve, S must be")
    print("the larger region: the punctured torus.")
    print("Therefore, the genus of S is g(S) = 1.")
    print("-" * 70)
    
    # Step 7: Final Calculation
    print("Step 7: We calculate the final genus g(M).")
    g_S = 1
    k = 1
    g_M = 2 * g_S + k - 1
    
    print(f"With g(S) = {g_S} and k = {k}, we plug these into the formula:")
    print(f"g(M) = 2 * g(S) + k - 1")
    print(f"The final calculation is:")
    print(f"{g_M} = 2 * {g_S} + {k} - 1")
    print(f"\nThe genus of the configuration space is {g_M}.")

solve_pentagon_genus()