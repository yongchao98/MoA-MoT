import math

def find_planar_distances():
    """
    Calculates the possible normalized distances (r <= 3) between hard spheres
    packed on a 2D plane, based on a hexagonal lattice structure.
    """
    # The limit for the normalized distance is r <= 3, so r^2 <= 9.
    limit_sq = 3**2

    # We use a set to store the unique squared distances to avoid duplicates.
    unique_distances_sq = set()

    # Iterate through a reasonable range of integer lattice coordinates (n, m).
    # If |n| > 3 or |m| > 3, r^2 will likely be > 9. A range of -3 to 3 is safe.
    for n in range(-4, 5):
        for m in range(-4, 5):
            # (0,0) is the central sphere itself, so we skip it.
            if n == 0 and m == 0:
                continue

            # For a hexagonal lattice, the squared distance is r^2 = n^2 + nm + m^2.
            r_sq = n**2 + n * m + m**2

            if r_sq <= limit_sq:
                unique_distances_sq.add(r_sq)

    # Sort the unique squared distances found.
    sorted_distances_sq = sorted(list(unique_distances_sq))

    print("The set of possible normalized distances r for r <= 3 are:")
    
    final_distances = []
    for r_sq in sorted_distances_sq:
        r = math.sqrt(r_sq)
        
        # Check if r is a whole number for cleaner output.
        # We use a small tolerance for floating point comparison.
        if abs(r - round(r)) < 1e-9:
            # It's an integer, print it directly.
            print(f"r = {r:.2f}")
        else:
            # It's not an integer, show the sqrt calculation.
            print(f"r = sqrt({r_sq}) = {r:.2f}")
        final_distances.append(f"{r:.2f}")
        
    # The final answer as a list of numbers as requested by the format <<<...>>>
    answer = ", ".join(final_distances)
    print(f"\n<<<{answer}>>>")


if __name__ == '__main__':
    find_planar_distances()