import math

def find_planar_distances():
    """
    Calculates the set of possible normalized distances (r) between hard spheres
    on a 2D plane for r <= 3.

    This is done by considering the two primary 2D lattice structures:
    hexagonal and square, as they represent the most probable ordered arrangements
    whose distances would appear as peaks in the radial distribution function g(r).
    """

    # Use a set to automatically store only unique distances.
    unique_distances = set()
    max_r = 3.0
    
    # We round distances to a high precision before adding to the set to handle
    # potential floating point inaccuracies and correctly identify unique values.
    precision = 5

    # A search range of 4 for lattice indices i and j is sufficient
    # to find all distances r <= 3.
    search_range = 4

    # --- Part 1: Hexagonal Lattice ---
    # The distance from the origin to a point (i,j) in a hexagonal lattice
    # (with unit distance between neighbors) is sqrt(i^2 + ij + j^2).
    for i in range(-search_range, search_range + 1):
        for j in range(-search_range, search_range + 1):
            if i == 0 and j == 0:
                continue
            
            # We check the squared distance first to avoid unnecessary sqrt() calls
            r_squared = float(i**2 + i * j + j**2)
            if r_squared <= max_r**2 + 1e-9: # Add tolerance for float comparison
                r = math.sqrt(r_squared)
                unique_distances.add(round(r, precision))

    # --- Part 2: Square Lattice ---
    # The distance from the origin to a point (i,j) in a square lattice
    # is sqrt(i^2 + j^2).
    for i in range(-search_range, search_range + 1):
        for j in range(-search_range, search_range + 1):
            if i == 0 and j == 0:
                continue

            r_squared = float(i**2 + j**2)
            if r_squared <= max_r**2 + 1e-9: # Add tolerance for float comparison
                r = math.sqrt(r_squared)
                unique_distances.add(round(r, precision))

    # --- Part 3: Format and Print Results ---
    # Convert the set to a sorted list for ordered output.
    sorted_distances = sorted(list(unique_distances))

    print("Based on the possible arrangements of hard spheres on a 2D plane (a hyperplane in 3D),")
    print("we consider the distances found in both hexagonal and square lattice packings.")
    print(f"The set of unique normalized distances r <= {max_r:.1f} is:")

    # Format the output as an equation showing each calculated distance.
    distance_strings = [f"{dist:.2f}" for dist in sorted_distances]
    print("r = {" + ", ".join(distance_strings) + "}")

if __name__ == '__main__':
    find_planar_distances()