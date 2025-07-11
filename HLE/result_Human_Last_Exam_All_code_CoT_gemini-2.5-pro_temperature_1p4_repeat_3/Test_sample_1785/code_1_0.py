def get_path_length(path):
    """Calculates the number of edges in a path."""
    # The number of edges is the number of vertices - 1 for a non-closed path.
    # For a closed path where the last vertex repeats the first, it's len(path) - 1.
    return len(path) - 1

def define_link_components():
    """
    Defines the vertex paths for a minimal nontrivial 3-component link.
    The link consists of a Hopf link (K1, K2) and a separate unknot (K3).
    """

    # K1: A 2x2 square in the z=0 plane.
    # Path follows the boundary of the square defined by (0,0,0) and (2,2,0).
    k1_path = [
        (0,0,0), (1,0,0), (2,0,0), (2,1,0), (2,2,0), (1,2,0), (0,2,0), (0,1,0), (0,0,0)
    ]

    # K2: A 2x2 square in the x=1 plane, positioned to link with K1.
    # Path follows the boundary of a square in its plane.
    k2_path = [
        (1,-1,-1), (1,0,-1), (1,1,-1), (1,1,0), (1,1,1), (1,0,1), (1,-1,1), (1,-1,0), (1,-1,-1)
    ]
    
    # K3: A 1x1 square (minimal unknot) placed far away to ensure disjoint vertices.
    k3_path = [
        (10,0,0), (11,0,0), (11,1,0), (10,1,0), (10,0,0)
    ]
    
    return k1_path, k2_path, k3_path

def main():
    """
    Constructs the minimal link and prints the total number of edges.
    """
    k1, k2, k3 = define_link_components()

    # Get the unique vertices for each component
    k1_vertices = set(k1)
    k2_vertices = set(k2)
    k3_vertices = set(k3)

    # Verify that the vertex sets are disjoint
    assert k1_vertices.isdisjoint(k2_vertices), "K1 and K2 share vertices"
    assert k1_vertices.isdisjoint(k3_vertices), "K1 and K3 share vertices"
    assert k2_vertices.isdisjoint(k3_vertices), "K2 and K3 share vertices"

    # Calculate the length of each knot (number of edges)
    len1 = get_path_length(k1)
    len2 = get_path_length(k2)
    len3 = get_path_length(k3)
    
    total_length = len1 + len2 + len3

    print("The minimum total number of edges is achieved with a Hopf link and a separate unknot.")
    print("Component 1 (K1) is a 2x2 loop with length: {}".format(len1))
    print("Component 2 (K2) is a 2x2 loop with length: {}".format(len2))
    print("Component 3 (K3) is a 1x1 loop with length: {}".format(len3))
    print("\nFinal Equation:")
    print("{} + {} + {} = {}".format(len1, len2, len3, total_length))

if __name__ == "__main__":
    main()
