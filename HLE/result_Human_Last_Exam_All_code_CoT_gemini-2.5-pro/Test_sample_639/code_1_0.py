import math

def calculate_distances():
    """
    Calculates the possible normalized distances (r <= 3) between two hard spheres
    on a 2D plane by considering triangular and square lattice arrangements.
    """
    # Use a set to store unique distances
    distances = set()
    
    # Define the limit
    r_max = 3.0
    r_max_squared = r_max * r_max

    # Maximum integer lattice steps to check
    # For r^2 = m^2 + n^2 <= 9, m,n <= 3
    # For r^2 = m^2 + mn + n^2 <= 9, m,n <= 3
    max_coord = 3 

    # 1. Triangular (hexagonal) lattice
    # r^2 = m^2 + mn + n^2
    for m in range(max_coord + 1):
        for n in range(max_coord + 1):
            # Skip the origin point
            if m == 0 and n == 0:
                continue
            
            r_squared = m*m + m*n + n*n
            if r_squared <= r_max_squared:
                distances.add(math.sqrt(r_squared))

    # 2. Square lattice
    # r^2 = m^2 + n^2
    for m in range(max_coord + 1):
        for n in range(max_coord + 1):
            # Skip the origin point
            if m == 0 and n == 0:
                continue
            
            r_squared = m*m + n*n
            if r_squared <= r_max_squared:
                distances.add(math.sqrt(r_squared))

    # Sort the unique distances
    sorted_distances = sorted(list(distances))

    # Print the results in the required format
    print("The set of normalized distances r <= 3 are:")
    for i, dist in enumerate(sorted_distances):
        print(f"r_{i+1} = {dist:.2f}")

if __name__ == "__main__":
    calculate_distances()