def solve_vector_puzzle():
    """
    Solves the vector puzzle by decrypting a message, rotating vectors,
    and summing them up.
    """
    message = "nggyunglydngraady"

    # Step 1 & 2: Decrypt message into a list of vectors based on a specific mapping
    mapping = {
        'g': (1, 0, 0),   # EAST
        'l': (-1, 0, 0),  # WEST
        'a': (0, 1, 0),   # NORTH
        'd': (0, -1, 0),  # SOUTH
        'r': (0, 0, 1),   # UP
        'y': (0, 0, -1)   # DOWN
    }

    initial_vectors = []
    for char in message:
        if char in mapping:
            initial_vectors.append(mapping[char])

    # Step 3: Rotate every 2nd vector (at odd-numbered positions, 2, 4, 6...)
    final_vectors = []
    for i, vec in enumerate(initial_vectors):
        # In a 1-based index, the 2nd vector is at index 1.
        # So we rotate for odd indices i (1, 3, 5...).
        if i % 2 != 0:
            x, y, z = vec
            rotated_vec = (x, z, -y)
            final_vectors.append(rotated_vec)
        else:
            final_vectors.append(vec)

    # Step 4: Sum all the vectors in the final list
    sum_x, sum_y, sum_z = 0, 0, 0
    for vec in final_vectors:
        sum_x += vec[0]
        sum_y += vec[1]
        sum_z += vec[2]

    final_sum = (sum_x, sum_y, sum_z)

    # Output the final equation showing all the vectors being summed
    equation_parts = [str(v) for v in final_vectors]
    equation_str = " + ".join(equation_parts)

    print("The final sum is derived from the following equation:")
    print(f"{equation_str} = {final_sum}")

solve_vector_puzzle()