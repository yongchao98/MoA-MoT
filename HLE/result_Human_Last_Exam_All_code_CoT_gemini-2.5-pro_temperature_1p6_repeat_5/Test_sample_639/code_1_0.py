import math

def calculate_and_print_distances():
    """
    Calculates and prints the possible normalized distances (r) between two hard spheres
    on a 2D plane for r <= 3.0, based on geometric arrangements.
    The distance r is normalized by the sphere diameter d.
    """
    print("Calculating possible distances 'r' for r <= 3...\n")

    # Use a set to automatically handle unique distances
    distances = set()

    # --- Distance 1: Direct Contact (r=1) ---
    # This is the most fundamental distance, where two spheres are touching.
    r1 = 1.0
    distances.add(r1)
    print(f"Arrangement: Direct contact of two spheres.")
    print(f"Calculation: By definition, the closest two non-overlapping spheres can be is one diameter.")
    print(f"r = 1.00\n")

    # --- Distance 2: Diagonal of a unit square (r=sqrt(2)) ---
    # This arrangement involves 4 spheres at the corners of a square.
    # The distance is between two diagonally opposite spheres.
    # d_ij^2 = d^2 + d^2 => (d_ij/d)^2 = 1^2 + 1^2 => r^2 = 2
    r2 = math.sqrt(2)
    distances.add(r2)
    print(f"Arrangement: Two spheres at opposite corners of a unit square of spheres.")
    print(f"Calculation: r = sqrt(1^2 + 1^2)")
    print(f"r = {r2:.2f}\n")

    # --- Distance 3: Second-nearest neighbor in hexagonal packing (r=sqrt(3)) ---
    # Two spheres (i, j) both touching a central sphere (k). The angle between the
    # lines from the center of k to i and j is 120 degrees in the next-neighbor case.
    # By the law of cosines: r^2 = 1^2 + 1^2 - 2*1*1*cos(120) = 2 - 2*(-0.5) = 3
    r3 = math.sqrt(3)
    distances.add(r3)
    print(f"Arrangement: Second-nearest neighbor in a hexagonal lattice.")
    print(f"Calculation: r = sqrt(3)")
    print(f"r = {r3:.2f}\n")

    # --- Distance 4: Collinear, separated by one sphere (r=2) ---
    # Three spheres in a straight line, i-k-j, all touching.
    # The distance is from the center of i to the center of j.
    r4 = 2.0
    distances.add(r4)
    print(f"Arrangement: Three spheres in a straight line, all touching.")
    print(f"Calculation: r = 1 + 1")
    print(f"r = {r4:.2f}\n")

    # --- Distance 5: From square lattice vector (2,1) (r=sqrt(5)) ---
    # Position vector is 2 units in x and 1 unit in y.
    # r^2 = 2^2 + 1^2 = 5
    r5 = math.sqrt(5)
    distances.add(r5)
    print(f"Arrangement: Relative position (2,1) in a square lattice.")
    print(f"Calculation: r = sqrt(2^2 + 1^2)")
    print(f"r = {r5:.2f}\n")

    # --- Distance 6: From hexagonal lattice (r=sqrt(7)) ---
    # Corresponds to a lattice vector of (2*e1 + 1*e2) in a hexagonal lattice
    # where e1=(1,0) and e2=(0.5, sqrt(3)/2). Vector is (2.5, sqrt(3)/2).
    # r^2 = 2.5^2 + (sqrt(3)/2)^2 = 6.25 + 0.75 = 7
    r6 = math.sqrt(7)
    distances.add(r6)
    print(f"Arrangement: Next-neighbor distance in a hexagonal lattice.")
    print(f"Calculation: r = sqrt(7)")
    print(f"r = {r6:.2f}\n")

    # --- Distance 7: From square lattice vector (2,2) (r=sqrt(8)) ---
    # Diagonal of a 2x2 square of spheres.
    # r^2 = 2^2 + 2^2 = 8
    r7 = math.sqrt(8)
    distances.add(r7)
    print(f"Arrangement: Diagonal of a 2x2 square of spheres.")
    print(f"Calculation: r = sqrt(2^2 + 2^2)")
    print(f"r = {r7:.2f}\n")
    
    # --- Distance 8: Collinear, separated by two spheres (r=3) ---
    # Four spheres in a straight line, i-k-l-j, all touching.
    r8 = 3.0
    distances.add(r8)
    print(f"Arrangement: Four spheres in a straight line, all touching.")
    print(f"Calculation: r = 1 + 1 + 1")
    print(f"r = {r8:.2f}\n")

    # --- Final Result ---
    # Sort the distances and print the final set.
    sorted_distances = sorted(list(distances))
    
    print("---------------------------------------------------------")
    print("The final set of possible distances for r <= 3 is:")
    final_distances_str = [f"{d:.2f}" for d in sorted_distances]
    print(final_distances_str)

calculate_and_print_distances()
<<<['1.00', '1.41', '1.73', '2.00', '2.24', '2.65', '2.83', '3.00']>>>