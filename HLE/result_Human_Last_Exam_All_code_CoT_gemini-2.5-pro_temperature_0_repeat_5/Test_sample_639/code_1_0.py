import math

def find_planar_distances():
    """
    Calculates the possible normalized distances (r <= 3) between hard spheres
    arranged in a 2D hexagonal lattice.
    """
    # The squared distance r^2 on a hexagonal lattice is given by n^2 + m^2 + nm.
    # We are looking for distances r <= 3, which means r^2 <= 9.
    # We can iterate through a range of integers n and m to find all possible distances.
    # A range of -3 to 3 for n and m is sufficient.
    
    max_r_squared = 9
    distances = set()

    for n in range(-3, 4):
        for m in range(-3, 4):
            # Skip the origin (distance to self)
            if n == 0 and m == 0:
                continue

            r_squared = n**2 + m**2 + n * m

            if r_squared <= max_r_squared:
                r = math.sqrt(r_squared)
                distances.add(r)

    # Sort the unique distances and print them
    sorted_distances = sorted(list(distances))
    
    print("The set of possible distances r for r <= 3 are:")
    for r in sorted_distances:
        print(f"{r:.2f}")

find_planar_distances()