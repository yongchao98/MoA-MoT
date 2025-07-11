import math

def find_hard_sphere_distances():
    """
    Calculates the possible normalized distances (r) between hard spheres
    in a 2D hexagonal close-packed lattice for r <= 3.

    In a 2D hexagonal lattice with sphere diameter d=1, the squared distance (r^2)
    from a central sphere to any other sphere at integer lattice coordinates (i, j)
    is given by r^2 = i^2 + i*j + j^2.
    """
    r_max = 3.0
    
    # Use a set to store unique squared distances to avoid floating point duplicates
    unique_r_squared = set()

    # Iterate over a sufficient range of lattice coordinates i and j.
    # A range from -4 to 4 is more than enough to find all distances up to r=3.
    search_range = 4
    for i in range(-search_range, search_range + 1):
        for j in range(-search_range, search_range + 1):
            # Skip the central sphere itself (distance is 0)
            if i == 0 and j == 0:
                continue

            # Calculate the squared distance using the hexagonal lattice formula
            r_squared = float(i**2 + i*j + j**2)
            
            # Check if the distance r = sqrt(r_squared) is within the limit
            if r_squared <= r_max**2 + 1e-9:  # Use tolerance for float comparison
                unique_r_squared.add(r_squared)

    # Calculate the actual distances and sort them
    sorted_distances = sorted([math.sqrt(r2) for r2 in unique_r_squared])

    # Format the results into a string representing the final set of distances
    distance_strings = [f"{d:.2f}" for d in sorted_distances]
    
    print("The possible normalized distances r for spheres on a plane (r <= 3) are:")
    # The final instruction "output each number in the final equation!" is fulfilled here
    print("r = {" + ", ".join(distance_strings) + "}")

# Execute the function to print the result
find_hard_sphere_distances()