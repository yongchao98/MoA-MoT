def solve_vector_puzzle():
    """
    Solves the vector puzzle by decrypting the message, generating vectors,
    rotating every second vector, and summing them up.
    """
    message = "nggyunglydngraady"

    # Define the mapping from direction characters to 3D vectors.
    # Z is vertical, so U/D are on the z-axis.
    # N/S are on the y-axis, and L/R are on the x-axis.
    direction_map = {
        'n': (0, 1, 0),
        'u': (0, 0, 1),
        'd': (0, 0, -1),
        'l': (-1, 0, 0),
        'r': (1, 0, 0)
    }

    # 1. Decrypt the message by filtering for valid direction characters
    # and create a list of corresponding vectors.
    vectors = []
    for char in message:
        if char in direction_map:
            vectors.append(direction_map[char])

    # Initialize the total sum vector.
    total_x, total_y, total_z = 0, 0, 0

    # 2. Iterate through the vectors, rotating every second one and summing them.
    for i, vec in enumerate(vectors):
        # The vector sequence is 1-based, so we check even positions
        # using the 0-based index `i`.
        # The 2nd, 4th, 6th... vectors correspond to indices 1, 3, 5...
        is_even_vector = (i + 1) % 2 == 0

        current_vec = vec
        if is_even_vector:
            # Rotate the vector 90 degrees clockwise around the x-axis.
            # The transformation is (x, y, z) -> (x, z, -y).
            x, y, z = vec
            current_vec = (x, z, -y)

        # 3. Add the current vector (rotated or not) to the total sum.
        total_x += current_vec[0]
        total_y += current_vec[1]
        total_z += current_vec[2]

    # 4. Print the final sum in the specified format.
    # The instruction "output each number in the final equation" is interpreted
    # as printing the final components of the vector.
    print(f"({total_x},{total_y},{total_z})")

solve_vector_puzzle()
<<<(0, 2, -2)>>>