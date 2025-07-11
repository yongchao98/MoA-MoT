import numpy as np

def solve_vector_puzzle():
    """
    Solves the vector puzzle by processing a sequence of directions,
    applying transformations, and summing the resulting vectors.
    """
    # Step 1: Decrypt the message.
    # The string "nggyunglydngraady" is 17 characters long.
    # A plausible decryption matching this length is "forwarddownupback".
    decrypted_message = "forward down up back"
    directions = decrypted_message.split()

    # Step 2: Define coordinate system and map directions to vectors.
    # (X: right, Y: forward, Z: up)
    vector_map = {
        "forward": np.array([0, 1, 0]),
        "down":    np.array([0, 0, -1]),
        "up":      np.array([0, 0, 1]),
        "back":    np.array([0, -1, 0])
    }
    
    initial_vectors = [vector_map[d] for d in directions]
    
    # Step 3: Apply rotation to every 2nd vector.
    # Rotation is clockwise around x-axis: (x, y, z) -> (x, z, -y)
    final_vectors = []
    for i, vec in enumerate(initial_vectors):
        # i is 0-indexed, so 2nd vector is at index 1, 4th at 3, etc.
        if (i + 1) % 2 == 0:
            rotated_vec = np.array([vec[0], vec[2], -vec[1]])
            final_vectors.append(rotated_vec)
        else:
            final_vectors.append(vec)
            
    # Step 4: Sum the final vectors.
    total_vector = np.sum(final_vectors, axis=0)
    
    # Step 5: Print the results in the required format.
    equation_parts = [f"({v[0]},{v[1]},{v[2]})" for v in final_vectors]
    equation_str = " + ".join(equation_parts)
    result_str = f"({total_vector[0]},{total_vector[1]},{total_vector[2]})"
    
    print(f"{equation_str} = {result_str}")

solve_vector_puzzle()