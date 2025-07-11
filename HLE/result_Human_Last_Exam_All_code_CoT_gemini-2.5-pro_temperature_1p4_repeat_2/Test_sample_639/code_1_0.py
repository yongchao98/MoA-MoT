import math

def find_packing_distances():
    """
    Calculates the set of possible normalized distances between hard spheres
    packed on a 2D plane for distances r <= 3.

    This is done by considering the two primary regular 2D packings:
    hexagonal and square lattices.
    """
    # Use a set to automatically store unique distances
    distances = set()

    # Define a search range for the integer lattice coordinates (n, m)
    # A range up to 5 is more than sufficient to find all distances r <= 3.
    search_limit = 5
    search_range = range(-search_limit, search_limit + 1)

    # Iterate through all combinations of lattice coordinates
    for n in search_range:
        for m in search_range:
            # Skip the origin particle at (0, 0)
            if n == 0 and m == 0:
                continue

            # --- Case 1: Hexagonal Lattice ---
            # The squared distance is r^2 = n^2 + nm + m^2
            r_squared_hex = n**2 + n * m + m**2
            # We are interested in r <= 3, which is equivalent to r^2 <= 9
            if r_squared_hex <= 9:
                r_hex = math.sqrt(r_squared_hex)
                distances.add(r_hex)

            # --- Case 2: Square Lattice ---
            # The squared distance is r^2 = n^2 + m^2
            r_squared_sq = n**2 + m**2
            # We are interested in r <= 3, which is equivalent to r^2 <= 9
            if r_squared_sq <= 9:
                r_sq = math.sqrt(r_squared_sq)
                distances.add(r_sq)

    # Convert the set to a sorted list for ordered output
    sorted_distances = sorted(list(distances))

    # Prepare the final output string
    output_parts = [
        "The set of possible distances r for r <= 3, normalized by the sphere diameter, are given by the equation:",
        "r = {" + ", ".join([f"{d:.2f}" for d in sorted_distances]) + "}"
    ]
    
    print("\n".join(output_parts))


if __name__ == "__main__":
    find_packing_distances()
    # The expected numerical values are the square roots of {1, 2, 3, 4, 5, 7, 8, 9}
    # which corresponds to {1.00, 1.41, 1.73, 2.00, 2.24, 2.65, 2.83, 3.00}
    final_numbers = [1.00, 1.41, 1.73, 2.00, 2.24, 2.65, 2.83, 3.00]
    final_answer = ",".join(map(str, final_numbers))
    # print(f'<<<{final_answer}>>>') # For internal check