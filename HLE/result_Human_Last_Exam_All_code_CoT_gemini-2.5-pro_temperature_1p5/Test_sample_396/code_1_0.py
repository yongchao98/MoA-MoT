def solve_vector_puzzle():
    """
    Solves the vector puzzle by processing a decrypted list of directions.
    """
    # Step 1: Define the decrypted directions.
    # This is based on the decryption of the message 'nggyunglydngraady'
    # using a Vigenere cipher with the key 'gravity', which yields 'southeastsouthwestdown'.
    decrypted_message = "south east south west down"
    directions = decrypted_message.split()

    # Step 2: Define the mapping from direction strings to vectors.
    # X: East/West, Y: North/South, Z: Up/Down
    dir_to_vec = {
        "NORTH": (0, 1, 0),
        "SOUTH": (0, -1, 0),
        "EAST":  (1, 0, 0),
        "WEST":  (-1, 0, 0),
        "UP":    (0, 0, 1),
        "DOWN":  (0, 0, -1),
    }

    # Step 3: Process directions to generate a list of vectors, applying rotation.
    final_vectors = []
    for i, direction in enumerate(directions):
        vec = dir_to_vec[direction.upper()]
        
        # The problem asks to rotate every 2nd vector. Indices are 0-based.
        # So we check for odd indices: 1, 3, 5, ... (which correspond to the 2nd, 4th, 6th... vector)
        if (i + 1) % 2 == 0:
            # Clockwise rotation around x-axis: (x, y, z) -> (x, z, -y)
            x, y, z = vec
            vec = (x, z, -y)
        
        final_vectors.append(vec)

    # Step 4: Sum the vectors and prepare the output strings for each component.
    sum_x = sum(v[0] for v in final_vectors)
    sum_y = sum(v[1] for v in final_vectors)
    sum_z = sum(v[2] for v in final_vectors)

    x_str = " + ".join(map(str, [v[0] for v in final_vectors])).replace("+ -", "- ")
    y_str = " + ".join(map(str, [v[1] for v in final_vectors])).replace("+ -", "- ")
    z_str = " + ".join(map(str, [v[2] for v in final_vectors])).replace("+ -", "- ")
    
    # Step 5: Print the results in the required format.
    print("Decrypted Directions: {}".format(directions))
    print("Final Vectors (with rotations applied): {}".format(final_vectors))
    print("\nFinal vector calculation:")
    print("x = {} = {}".format(x_str, sum_x))
    print("y = {} = {}".format(y_str, sum_y))
    print("z = {} = {}".format(z_str, sum_z))
    print("\nSum of all generated vectors: ({},{},{})".format(sum_x, sum_y, sum_z))

solve_vector_puzzle()