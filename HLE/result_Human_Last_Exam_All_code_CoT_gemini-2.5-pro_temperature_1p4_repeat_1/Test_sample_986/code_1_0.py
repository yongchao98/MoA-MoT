import itertools

def is_adjacent(v1, v2):
    """
    Checks if two vertices v1 and v2 in the graph X are adjacent.
    A vertex is a tuple (x, y) with x <= y.
    Adjacency rule: y1 = x2 or y2 = x1.
    """
    x1, y1 = v1
    x2, y2 = v2
    return y1 == x2 or y2 == x1

def is_clique(vertex_set):
    """
    Checks if a given set of vertices forms a clique.
    A clique is a set where every two distinct vertices are adjacent.
    """
    for v1, v2 in itertools.combinations(vertex_set, 2):
        if not is_adjacent(v1, v2):
            return False
    return True

def solve():
    """
    This function finds and demonstrates the clique number of the graph X.
    Based on mathematical analysis, the clique number is 3.
    We demonstrate this by providing an example of a 3-clique.
    """
    # A set of 3 vertices in the graph X.
    # Each vertex (x, y) must satisfy x <= y.
    # (-1, 0) satisfies -1 <= 0.
    # (0, 0) satisfies 0 <= 0.
    # (0, 1) satisfies 0 <= 1.
    example_clique = {(-1, 0), (0, 0), (0, 1)}

    # Verify that this set is indeed a clique.
    if is_clique(example_clique):
        clique_number = len(example_clique)
        print(f"An example of a maximum clique is the set of vertices: {example_clique}")
        # The prompt requires outputting each number in the final equation.
        # We interpret this as showing the calculation of the size of the clique.
        print(f"The clique consists of the following {clique_number} vertices:")
        for v in example_clique:
            print(f"  Vertex: {v}")
        
        calculation_string = " + ".join(["1"] * clique_number)
        print(f"The size of this clique is {calculation_string} = {clique_number}.")
        print(f"The clique number of X is {clique_number}.")
    else:
        # This part should not be reached based on our analysis.
        print("The provided set is not a clique.")

solve()