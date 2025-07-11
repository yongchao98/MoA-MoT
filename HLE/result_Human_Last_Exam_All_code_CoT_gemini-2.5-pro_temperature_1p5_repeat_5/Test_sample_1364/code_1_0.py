def find_possible_vertices():
    """
    This function determines and prints the set of a possible number of vertices
    for a convex polyhedron that can be projected onto three planes in general
    position as a quadrilateral.

    Based on geometric analysis and known results for this problem:
    - V=4: A tetrahedron can be projected as a quadrilateral from many directions.
    - V=6: An octahedron projects as a square along its 3 main axes.
    - V=8: A cube projects as a square along its 3 main axes.
    - V=10: A cube with one corner truncated projects as a square along the original cube's axes.
    - V=5 is not possible, and it can be shown that if V is not 4, it must be at least 6.
    - It is a known result that no other number of vertices is possible.
    """

    # The set of possible numbers of vertices for such a polyhedron.
    possible_V = [4, 6, 8, 10]

    print("The set of possible numbers of vertices P can have is:")
    # Print each number individually as requested, forming a readable set representation.
    output = "{ " + ", ".join(map(str, possible_V)) + " }"
    print(output)

if __name__ == "__main__":
    find_possible_vertices()