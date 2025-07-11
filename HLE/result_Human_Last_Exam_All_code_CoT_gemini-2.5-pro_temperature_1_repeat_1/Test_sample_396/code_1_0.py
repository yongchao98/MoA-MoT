import collections

def solve_vector_sum():
    """
    Calculates the sum of vectors based on a decrypted message,
    applying specific rotation rules.
    """
    # Step 1: Assume the decrypted sequence of movements.
    # This is based on the reasoning that the provided ciphertext is ambiguous.
    directions = ["south", "south", "west", "up"]

    # Step 2: Define the vector for each direction.
    # X=East, Y=North, Z=Up
    vector_map = {
        "north": (0, 1, 0),
        "south": (0, -1, 0),
        "east":  (1, 0, 0),
        "west":  (-1, 0, 0),
        "up":    (0, 0, 1),
        "down":  (0, 0, -1),
    }

    final_vectors = []
    
    # Step 3: Process each direction, applying rotation to every 2nd vector.
    for i, direction in enumerate(directions):
        vector = collections.namedtuple('Vector', ['x', 'y', 'z'])(*vector_map[direction])
        
        # Vectors are 1-indexed for the problem, so check (i+1) for evenness.
        is_even_vector = (i + 1) % 2 == 0

        if is_even_vector:
            # Apply clockwise rotation around the x-axis: (x, y, z) -> (x, z, -y)
            rotated_vector = vector_map[direction]
            final_vector = (rotated_vector[0], rotated_vector[2], -rotated_vector[1])
            final_vectors.append(final_vector)
        else:
            final_vectors.append(vector_map[direction])

    # Step 4: Calculate the sum of all generated vectors.
    sum_x = sum(v[0] for v in final_vectors)
    sum_y = sum(v[1] for v in final_vectors)
    sum_z = sum(v[2] for v in final_vectors)

    # Step 5: Format the output to show each number in the final equation.
    x_components = [str(v[0]) for v in final_vectors]
    y_components = [str(v[1]) for v in final_vectors]
    z_components = [str(v[2]) for v in final_vectors]
    
    # To make the plus signs look good for negative numbers
    x_str = "+".join(x_components).replace("+-", "-")
    y_str = "+".join(y_components).replace("+-", "-")
    z_str = "+".join(z_components).replace("+-", "-")

    print(f"({x_str}, {y_str}, {z_str}) = ({sum_x},{sum_y},{sum_z})")

solve_vector_sum()
<<<(-1,0,1)>>>