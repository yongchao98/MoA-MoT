def find_minimum_vertices():
    """
    Calculates and explains the minimum number of vertices in a two-loop
    Feynman diagram for a phi-4 interacting scalar field theory.
    """
    # Number of loops
    L = 2

    print("To find the minimum number of vertices in a two-loop Feynman diagram,")
    print("we use the general relations for a scalar phi-4 theory.\n")

    print(f"First, we have the topological formula relating loops (L), internal lines (I), and vertices (V):")
    print(f"L = I - V + 1")
    print(f"For L = {L} loops, this becomes: {L} = I - V + 1, which rearranges to I = V + 1.\n")

    print("Second, for a phi-4 theory, each vertex connects 4 lines. This gives the relation:")
    print("4 * V = 2 * I + E, where E is the number of external lines.\n")

    print("Now, we substitute the first equation into the second:")
    print("4 * V = 2 * (V + 1) + E")
    print("4 * V = 2 * V + 2 + E\n")

    print("Solving this equation for V gives us:")
    print("2 * V = 2 + E")
    print("V = (2 + E) / 2\n")

    print("To find the minimum number of vertices (V), we must use the minimum possible")
    print("number of external lines (E). A valid Feynman diagram can have zero external")
    print("lines (E=0). These are known as 'vacuum diagrams'.\n")

    # Minimum number of external lines
    E_min = 0

    # Calculate the minimum number of vertices
    V_min = (2 + E_min) / 2

    print("Using the minimum value E = 0, the minimum number of vertices is:")
    print(f"V_min = (2 + {E_min}) / 2 = {int(V_min)}")
    print("\nThis diagram is the 'figure-eight' vacuum diagram, which has 1 vertex, 2 internal lines, and forms 2 loops.")

find_minimum_vertices()