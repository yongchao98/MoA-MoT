import math

def calculate_hard_sphere_distances():
    """
    Calculates the possible normalized distances (r) between two hard spheres
    on a 2D plane for r <= 3, considering hexagonal and square packing configurations.
    """
    # Use a set to store unique distances to avoid duplicates
    unique_distances = set()
    
    # Maximum integer coordinate to check for lattice points (n, m)
    # Since 3*3 = 9, a max coordinate of 3 is sufficient.
    max_coord = 3
    r_max = 3.0

    print("Finding possible normalized distances r <= 3.00 for 2D packed spheres.")
    print("-" * 60)

    # --- 1. Hexagonal Lattice ---
    # The squared distance is r^2 = n^2 + nm + m^2 for lattice vectors (n, m)
    for n in range(max_coord + 1):
        for m in range(max_coord + 1):
            if n == 0 and m == 0:
                continue  # Skip the reference sphere itself
            
            r_squared = n**2 + n*m + m**2
            r = math.sqrt(r_squared)
            
            if r <= r_max:
                unique_distances.add(r)

    # --- 2. Square Lattice ---
    # The squared distance is r^2 = n^2 + m^2 for lattice vectors (n, m)
    for n in range(max_coord + 1):
        for m in range(max_coord + 1):
            if n == 0 and m == 0:
                continue # Skip the reference sphere itself

            r_squared = n**2 + m**2
            r = math.sqrt(r_squared)
            
            if r <= r_max:
                unique_distances.add(r)

    # Sort the final list of unique distances
    sorted_distances = sorted(list(unique_distances))
    
    print("The set of unique distances found is:")
    final_list_str = []
    for r in sorted_distances:
        # Format to two decimal places and print
        formatted_r = f"{r:.2f}"
        print(formatted_r)
        final_list_str.append(formatted_r)

    # Required final answer format
    final_answer = ", ".join(final_list_str)
    print("\nFinal Answer for Submission:")
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    calculate_hard_sphere_distances()
