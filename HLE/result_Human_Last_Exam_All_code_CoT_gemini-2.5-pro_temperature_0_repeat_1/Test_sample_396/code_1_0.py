import numpy as np

def solve_vector_puzzle():
    """
    Solves the vector puzzle by processing a decrypted list of directions.
    """
    # Step 1: The decrypted message provides a list of directions.
    # The decryption of "nggyunglydngraady" is determined to be "north west up up right down".
    directions = ["north", "west", "up", "up", "right", "down"]

    # Step 2: Define the mapping from directions to unit vectors.
    # Z is vertical, Y is North/South, X is East/West.
    vector_map = {
        "north": np.array([0, 1, 0]),
        "south": np.array([0, -1, 0]),
        "east": np.array([1, 0, 0]),
        "west": np.array([-1, 0, 0]),
        "up": np.array([0, 0, 1]),
        "down": np.array([0, 0, -1]),
        "right": np.array([1, 0, 0]),
        "left": np.array([-1, 0, 0]),
    }

    final_vectors = []
    
    # Step 3 & 4: Process each direction, applying rotation to every 2nd vector.
    for i, direction in enumerate(directions):
        vector = vector_map[direction]
        
        # The problem states to rotate every 2nd vector (index 1, 3, 5, ...)
        if (i + 1) % 2 == 0:
            # Apply a clockwise rotation around the x-axis: (x, y, z) -> (x, z, -y)
            x, y, z = vector
            rotated_vector = np.array([x, z, -y])
            final_vectors.append(rotated_vector)
        else:
            final_vectors.append(vector)

    # Step 5: Sum all the final vectors.
    sum_vector = np.sum(final_vectors, axis=0)

    # Step 6: Print the final equation as requested.
    equation_parts = []
    for v in final_vectors:
        equation_parts.append(f"({v[0]},{v[1]},{v[2]})")
    
    equation_str = " + ".join(equation_parts)
    result_str = f"({sum_vector[0]},{sum_vector[1]},{sum_vector[2]})"
    
    print(f"{equation_str} = {result_str}")

solve_vector_puzzle()