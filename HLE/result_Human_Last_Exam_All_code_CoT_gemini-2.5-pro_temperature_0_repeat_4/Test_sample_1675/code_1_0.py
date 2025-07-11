def solve_max_points_problem():
    """
    This script explains the step-by-step solution to find the maximum
    number of points n.
    """

    print("Let R, G, Y be the number of red, green, and yellow points, respectively.")
    print("The total number of points is n = R + G + Y.")
    print("\nThe three conditions are:")
    print("1. Any triangle of 3 red points contains a green point.")
    print("2. Any triangle of 3 green points contains a yellow point.")
    print("3. Any triangle of 3 yellow points contains a red point.")
    print("-" * 30)

    print("\nStep 1: Rephrasing the conditions using convex hulls.")
    print("Let P_R, P_G, P_Y be the sets of red, green, and yellow points.")
    print("Let's analyze condition 1. If R >= 3, consider any red point r. If r is outside the convex hull of all other points (P_G U P_Y), we could form a red triangle containing r that is separated from all green points, which is a contradiction.")
    print("Therefore, if R >= 3, every red point must be inside the convex hull of the green and yellow points. So, P_R is a subset of CH(P_G U P_Y).")
    print("By symmetry, the same logic applies to all colors:")
    print("- If R >= 3, then P_R is a subset of CH(P_G U P_Y).")
    print("- If G >= 3, then P_G is a subset of CH(P_R U P_Y).")
    print("- If Y >= 3, then P_Y is a subset of CH(P_R U P_G).")
    print("-" * 30)

    print("\nStep 2: Proving that not all color sets can have >= 3 points.")
    print("Let H be the convex hull of all n points, H = CH(P_R U P_G U P_Y).")
    print("The vertices of H must be points from P_R, P_G, or P_Y.")
    print("Suppose R >= 3. Then P_R is a subset of CH(P_G U P_Y). This implies that no red point can be a vertex of H. (A vertex of a convex hull cannot be a convex combination of other points in the set).")
    print("Similarly, if G >= 3, no green point can be a vertex of H.")
    print("And if Y >= 3, no yellow point can be a vertex of H.")
    print("If we assume R>=3, G>=3, and Y>=3 all hold, then no point of any color can be a vertex of H. This is a contradiction, as H must have vertices.")
    print("Therefore, at least one color set must have fewer than 3 points.")
    print("-" * 30)

    print("\nStep 3: Proving that at least two color sets must have <= 2 points.")
    print("From Step 2, at least one color count is at most 2. Let's assume Y <= 2.")
    print("Now, suppose both R >= 3 and G >= 3.")
    print("As shown before, this means no red or green point can be a vertex of the overall convex hull H.")
    print("So, all vertices of H must be yellow. The number of vertices of H must be at least 3 (since no three points are collinear).")
    print("But we assumed Y <= 2, so there are at most 2 yellow points. This means H can have at most 2 vertices, which is a contradiction.")
    print("So, the assumption that R >= 3 and G >= 3 must be false.")
    print("This means that if Y <= 2, then either R <= 2 or G <= 2.")
    print("Conclusion: At least two of the three color sets must have 2 or fewer points.")
    print("-" * 30)

    print("\nStep 4: Maximizing the size of the third color set.")
    print("Let's assume G <= 2 and Y <= 2. To maximize n, we should take G=2 and Y=2.")
    print("The conditions on green and yellow triangles are now vacuously true.")
    print("We only need to satisfy Condition 1: Any triangle of 3 red points must contain a green point.")
    print("Let the two green points be g1 and g2. Consider the line L passing through g1 and g2.")
    print("The line L divides the plane into two open half-planes, H1 and H2.")
    print("If we choose 3 red points from H1, the triangle they form will lie entirely in H1 and cannot contain g1 or g2. This is a contradiction.")
    print("Therefore, there can be at most 2 red points in H1.")
    print("Similarly, there can be at most 2 red points in H2.")
    print("(No red point can lie on L because of the 'no three points collinear' rule).")
    print("So, the maximum number of red points is R_max = 2 + 2 = 4.")
    print("-" * 30)

    print("\nStep 5: Calculating the maximum value of n.")
    print("We found the maximum number of points is achieved with a distribution of (R, G, Y) = (4, 2, 2) or its permutations.")
    R, G, Y = 4, 2, 2
    n = R + G + Y
    print(f"The maximum value of n is R + G + Y = {R} + {G} + {Y} = {n}.")
    print("\nA possible construction for (4, 2, 2):")
    print("Place g1 at (0, 1) and g2 at (0, -1).")
    print("Place two red points in the right half-plane, e.g., r1=(1,1), r2=(1,0).")
    print("Place two red points in the left half-plane, e.g., r3=(-1,1), r4=(-1,0).")
    print("Place the two yellow points anywhere else, e.g., y1=(10,10), y2=(-10,-10).")
    print("Any red triangle must pick vertices from both half-planes, and will thus contain the green points on the y-axis.")

if __name__ == '__main__':
    solve_max_points_problem()
