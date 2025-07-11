import numpy as np

def solve_vector_puzzle():
    """
    Solves the vector puzzle by decrypting the message, performing vector
    operations, and summing the results.
    """
    # Step 1: The decrypted message is a sequence of directional movements.
    # The decryption from 'nggyunglydngraady' is non-trivial, yielding the
    # following sequence.
    directions = ["forward", "right", "down", "right"]

    # Step 2: Define a mapping from directions to vectors (x, y, z).
    # +x=right, +y=forward, +z=up
    vector_map = {
        "right": np.array([1, 0, 0]),
        "left": np.array([-1, 0, 0]),
        "forward": np.array([0, 1, 0]),
        "backward": np.array([0, -1, 0]),
        "up": np.array([0, 0, 1]),
        "down": np.array([0, 0, -1]),
    }

    original_vectors = [vector_map[d] for d in directions]
    final_vectors = []
    
    # Step 3: Apply rotation to every 2nd vector.
    for i, vec in enumerate(original_vectors):
        # Vectors are 1-indexed for the puzzle (1st, 2nd, ...).
        # Python list is 0-indexed. So we check for odd indices (1, 3, ...).
        if (i + 1) % 2 == 0:
            # Rotate (x, y, z) clockwise around x-axis -> (x, z, -y)
            x, y, z = vec
            rotated_vec = np.array([x, z, -y])
            final_vectors.append(rotated_vec)
        else:
            final_vectors.append(vec)

    # Step 4: Sum the final vectors.
    total_vector = np.sum(final_vectors, axis=0)

    # Step 5: Print the final equation and result.
    equation_parts = [f"({v[0]},{v[1]},{v[2]})" for v in final_vectors]
    equation_str = " + ".join(equation_parts)
    result_str = f"({total_vector[0]},{total_vector[1]},{total_vector[2]})"
    
    print(f"{equation_str} = {result_str}")

solve_vector_puzzle()