import numpy as np

def solve_vector_puzzle():
    """
    Solves the vector puzzle by decrypting the message, calculating vectors,
    applying rotations, and summing them up.
    """
    # Step 1 & 3: Decrypt message and create initial vectors.
    # The message 'nggyunglydngraady' is assumed to decrypt to this sequence.
    directions = ["right", "up", "left", "down", "forward", "backward"]
    
    # Mapping directions to vectors based on Z-vertical and unit movement
    vector_map = {
        "right":    np.array([1, 0, 0]),
        "left":     np.array([-1, 0, 0]),
        "forward":  np.array([0, 1, 0]),
        "backward": np.array([0, -1, 0]),
        "up":       np.array([0, 0, 1]),
        "down":     np.array([0, 0, -1]),
    }

    initial_vectors = [vector_map[d] for d in directions]
    
    final_vectors = []
    
    # Step 4: Rotate every 2nd vector
    for i, vec in enumerate(initial_vectors):
        # The problem asks to rotate the 2nd, 4th, etc., vectors (1-indexed)
        if (i + 1) % 2 == 0:
            x, y, z = vec
            # Clockwise rotation around x-axis: (x, y, z) -> (x, z, -y)
            rotated_vec = np.array([x, z, -y])
            final_vectors.append(rotated_vec)
        else:
            final_vectors.append(vec)

    # Step 5: Sum the final vectors
    sum_vector = np.sum(final_vectors, axis=0)
    
    # Step 6: Format and print the final equation
    vector_strings = [f"({v[0]},{v[1]},{v[2]})" for v in final_vectors]
    equation_str = " + ".join(vector_strings)
    result_str = f"({sum_vector[0]},{sum_vector[1]},{sum_vector[2]})"
    
    print(f"{equation_str} = {result_str}")

solve_vector_puzzle()