import math

def calculate_particle_distances():
    """
    Calculates the possible distances between hard spheres packed on a 2D plane
    for both square and hexagonal lattice arrangements, up to a normalized distance of 3.
    """
    # The maximum normalized distance to consider
    r_max = 3
    
    # Use a set to store unique distances to avoid duplicates
    unique_distances = set()

    # Iterate through a grid of integer indices (i, j).
    # A range up to r_max is sufficient, as i or j > r_max would result in r > r_max.
    for i in range(r_max + 1):
        for j in range(r_max + 1):
            # Skip the origin (distance to itself)
            if i == 0 and j == 0:
                continue

            # --- Case 1: Square Lattice ---
            # The squared distance is r^2 = i^2 + j^2
            r_sq_square = i**2 + j**2
            if r_sq_square <= r_max**2:
                distance_square = math.sqrt(r_sq_square)
                unique_distances.add(distance_square)

            # --- Case 2: Hexagonal Lattice ---
            # The squared distance is r^2 = i^2 + j^2 + i*j
            r_sq_hex = i**2 + j**2 + i * j
            if r_sq_hex <= r_max**2:
                distance_hex = math.sqrt(r_sq_hex)
                unique_distances.add(distance_hex)
    
    # Sort the collected distances in ascending order
    sorted_distances = sorted(list(unique_distances))

    print("The set of possible normalized distances (r) for r <= 3 are:")
    # Print each distance formatted to two decimal places
    for r in sorted_distances:
        print(f"r = {r:.2f}")

# Execute the function
calculate_particle_distances()

# The final answer as a list of numbers
final_answer = sorted(list({math.sqrt(i**2+j**2) for i in range(4) for j in range(4) if 0 < i**2+j**2 <= 9} | {math.sqrt(i**2+j**2+i*j) for i in range(4) for j in range(4) if 0 < i**2+j**2+i*j <= 9}))
# <<<1.00, 1.41, 1.73, 2.00, 2.24, 2.65, 2.83, 3.00>>>