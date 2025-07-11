import math

def find_coplanar_distances():
    """
    Calculates the set of possible normalized distances (r) between two hard spheres
    in a coplanar arrangement for r <= 3.

    This is done by considering the geometry of the two primary 2D crystal structures:
    hexagonal and square lattices.
    """
    # Use a set to automatically store only unique distances
    distances = set()

    # Define a search range for lattice indices. A range from -4 to 4 is sufficient
    # to find all distances where r <= 3.
    search_range = range(-4, 5)

    # --- Part 1: Hexagonal Lattice Distances ---
    # The squared distance r^2 in a hexagonal lattice is given by m^2 + mn + n^2,
    # where m and n are integer lattice indices.
    for m in search_range:
        for n in search_range:
            # Skip the origin particle (distance to itself is 0)
            if m == 0 and n == 0:
                continue
            
            r_squared = float(m**2 + m * n + n**2)
            r = math.sqrt(r_squared)
            
            # Check if the distance is within the desired range (r <= 3).
            # An epsilon (0.0001) is used for robust float comparison.
            if r <= 3.0001:
                distances.add(r)

    # --- Part 2: Square Lattice Distances ---
    # The squared distance r^2 in a square lattice is given by i^2 + j^2,
    # where i and j are integer lattice indices.
    for i in search_range:
        for j in search_range:
            # Skip the origin particle
            if i == 0 and j == 0:
                continue
                
            r_squared = float(i**2 + j**2)
            r = math.sqrt(r_squared)
            
            # Check if the distance is within the desired range
            if r <= 3.0001:
                distances.add(r)
                
    # --- Part 3: Process and Print Results ---
    # Convert the set to a list and sort it in ascending order
    sorted_distances = sorted(list(distances))
    
    print("The set of possible normalized distances (r) for r <= 3 is:")
    
    final_answer_list = []
    # Print each resulting distance formatted to two decimal places
    for r in sorted_distances:
        formatted_r = f"{r:.2f}"
        print(formatted_r)
        final_answer_list.append(formatted_r)
        
    # Format the final answer as a single string for the final output
    final_answer_string = ", ".join(final_answer_list)
    print(f"<<<{final_answer_string}>>>")

if __name__ == "__main__":
    find_coplanar_distances()