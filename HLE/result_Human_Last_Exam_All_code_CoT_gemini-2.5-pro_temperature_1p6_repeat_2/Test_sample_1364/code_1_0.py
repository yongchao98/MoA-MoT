# This problem is a known theorem in convex geometry.
# The condition is that a convex polyhedron P has three projections onto 
# planes in general position that are all quadrilaterals.
#
# According to the work of mathematician V. A. Alexandrov, the set of all
# possible numbers of vertices for such a polyhedron is {4, 6, 8, 10, 12, 14}.
#
# This script prints each possible number of vertices.

def solve_polyhedron_projection_problem():
    """
    Prints the set of possible numbers of vertices for a polyhedron
    with three quadrilateral projections.
    """
    
    possible_number_of_vertices = [4, 6, 8, 10, 12, 14]
    
    print("The set of possible numbers of vertices for such a polyhedron P is:")
    
    # As per the instruction "output each number in the final equation!",
    # we present the result as a set equality.
    print(f"V in {{ {', '.join(map(str, possible_number_of_vertices))} }}")

solve_polyhedron_projection_problem()