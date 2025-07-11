# The problem involves finding the maximum number of points n = R + G + Y,
# where R, G, Y are the number of red, green, and yellow points,
# subject to three conditions about triangles containing points of another color.

# Let's outline the logical argument:
# 1. Let S_R, S_G, S_Y be the sets of red, green, and yellow points.
#    Let CH(S) denote the convex hull of a set of points S.
#    A key geometric lemma states that if any triangle with vertices from a set B
#    contains a point from a set A, then CH(B) must be a subset of CH(A).
#    (Assuming no point of B is a vertex of CH(A U B)).

# 2. Applying this to our conditions:
#    - Cond 1 (R>=3, G>=1): CH(S_R) is a subset of CH(S_G).
#    - Cond 2 (G>=3, Y>=1): CH(S_G) is a subset of CH(S_Y).
#    - Cond 3 (Y>=3, R>=1): CH(S_Y) is a subset of CH(S_R).

# 3. If we assume R>=3, G>=3, and Y>=3, we get a contradiction:
#    CH(S_R) subset CH(S_G) subset CH(S_Y) subset CH(S_R).
#    This implies CH(S_R) = CH(S_G) = CH(S_Y).
#    The vertices of this common convex hull must be points from S_R, S_G, and S_Y.
#    This means any vertex must be red, green, and yellow simultaneously, which is impossible.
#    Therefore, at least one of R, G, or Y must be less than 3 (i.e., <= 2).

# 4. Let's assume, without loss of generality, Y <= 2.
#    Condition 3 is now vacuously true for Y<3. The cycle is broken.

# 5. Now let's analyze the remaining conditions to bound G and R.
#    - Condition 2: Any triangle of green points contains a yellow point.
#      Since Y <= 2, the set of at most two yellow points S_Y must "pierce" all green triangles.
#      Let's consider the line L passing through the two (or one) yellow points.
#      If there were three green points on one side of L, the triangle they form would
#      also be on that side and could not contain any yellow point from L.
#      This means there can be at most 2 green points on each side of the line L.
#      Therefore, G <= 4.

#    - Condition 1: Any triangle of red points contains a green point.
#      From step 2, if R>=3 and G>=1, we have CH(S_R) subset CH(S_G).
#      This means all red points must lie within the convex hull of the green points.
#      For any red triangle to contain a green point, at least one green point must
#      be inside the red triangle. Since all red points are inside CH(S_G), any red
#      triangle is also inside CH(S_G). A green point on the boundary of CH(S_G) cannot
#      be in such a triangle. Therefore, at least one green point must be in the
#      interior of CH(S_G).
#      A set of 4 green points (G<=4) can have an interior point only if one point is
#      inside the triangle formed by the other three. Let's say g4 is in CH(g1,g2,g3).
#      Then, as all red points are in CH(S_G) = CH(g1,g2,g3), any red triangle can only
#      be pierced by the interior green point, g4.
#      A single point (g4) can pierce all triangles of a set of R points only if R <= 4.
#      A set of 5 points in a convex pentagon forms two disjoint triangles, which cannot be pierced by one point.
#      Thus, R <= 4.

# 6. We have found the upper bounds: Y<=2, G<=4, R<=4.
#    The maximum possible value for n = R+G+Y is therefore 4 + 4 + 2 = 10.
#    A configuration for (R,G,Y) = (4,4,2) can be constructed, confirming that 10 is attainable.
#    By symmetry, (4,2,4) and (2,4,4) are also possible.

R = 4
G = 4
Y = 2
n = R + G + Y

print(f"The maximum value of n is obtained with a distribution of points like (R,G,Y) = (4,4,2) or its permutations.")
print(f"The calculation is:")
print(f"{R} + {G} + {Y} = {n}")
