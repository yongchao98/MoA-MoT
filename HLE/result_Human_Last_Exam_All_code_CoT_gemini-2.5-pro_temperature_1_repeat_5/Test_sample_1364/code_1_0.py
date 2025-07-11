def solve_polyhedron_problem():
    """
    This function prints the step-by-step reasoning to solve the problem
    about the possible number of vertices for a specific convex polyhedron.
    """
    
    explanation = """
# Step 1: Analyze the problem statement.
# The condition is that a 3D convex polyhedron P has three projections (onto three
# planes in general position) that are all quadrilaterals.
# This means for three linearly independent projection directions, the "shadow boundary"
# on the polyhedron is a 4-cycle of vertices.

# Step 2: Find the minimum possible number of vertices (V).

print("--- Finding the Lower Bound ---")

# V=4 (Tetrahedron)
print("V=4: A tetrahedron's projection is always a triangle. Not possible.")

# V=5 (e.g., Square Pyramid, Triangular Bipyramid)
print("V=5: A square pyramid has at most one projection that is a quadrilateral.")
print("       A triangular bipyramid can have quadrilateral projections, but not for three linearly independent directions.")
print("       Therefore, V=5 is not possible.")

# V=6 (Octahedron)
print("V=6: An octahedron with vertices at (+-1, 0, 0), (0, +-1, 0), (0, 0, +-1) works.")
print("       The projections onto the xy, yz, and xz planes are all squares (quadrilaterals).")
print("       The projection directions are orthogonal and thus linearly independent.")
print("       Therefore, V=6 is the minimum possible number of vertices.")

# Step 3: Show that all integers V >= 6 are possible using an inductive construction.

print("\\n--- Proving for all V >= 6 by Induction ---")

# Base Case:
print("Base Case: V=6 is possible, as shown with the octahedron.")

# Inductive Step:
print("Inductive Step: Assume a polyhedron P_k with k vertices (k>=6) satisfies the condition.")
print("1. For any three linearly independent directions, there must be a face F on P_k that is visible from all three directions.")
print("2. We can form a new polyhedron P_{k+1} by 'capping' this face F with a new vertex, v_new.")
print("3. This new polyhedron has k+1 vertices.")
print("4. Since F was visible from all three directions, v_new projects inside the existing quadrilateral projections.")
print("5. Thus, P_{k+1} also satisfies the condition.")

# Conclusion of proof:
print("Conclusion: Since V=6 works and we can always construct a valid polyhedron with V=k+1 from one with V=k, all integers V >= 6 are possible.")

# Final Answer
print("\\n---------------------------------")
print("The set of possible numbers of vertices is all integers greater than or equal to 6.")
print("Set = {6, 7, 8, 9, ...}")
print("---------------------------------")

"""
    # Execute the explanation printout
    exec(explanation)

solve_polyhedron_problem()