import math

def find_planar_distances():
    """
    Calculates the possible normalized distances (r <= 3) between hard spheres
    packed on a 2D plane.

    The calculation assumes the spheres form a hexagonal lattice, which is the densest
    packing in 2D. The squared distance r^2 from a lattice point at the origin
    to another lattice point (u, v) is given by r^2 = u^2 + u*v + v^2.
    """
    # We are looking for distances r <= 3, so r^2 <= 9.
    # A search range for u and v up to 3 is sufficient.
    search_range = 4
    
    # Use a dictionary to store one example (u, v) pair for each unique r^2 found.
    # This avoids storing duplicate distances and helps with the explanation.
    # Key: r_squared, Value: (u, v)
    distance_map = {}

    for u in range(search_range):
        for v in range(-search_range, search_range):
            # The case u=0, v=0 is the distance to the sphere itself (r=0), which we exclude.
            if u == 0 and v == 0:
                continue

            r_squared = float(u**2 + u*v + v**2)

            # We only care about distances where r <= 3, so r^2 <= 9.
            if r_squared <= 9:
                # If we haven't seen this distance before, store it with its (u,v) generators.
                if r_squared not in distance_map:
                    distance_map[r_squared] = (u, v)

    # Sort the unique squared distances in ascending order.
    sorted_r_squared = sorted(distance_map.keys())

    print("The possible normalized distances r for r <= 3 are:\n")
    
    final_distances = []
    # Print the calculation for each unique distance found.
    for r_sq in sorted_r_squared:
        u, v = distance_map[r_sq]
        r = math.sqrt(r_sq)
        
        # Display the formula and the result for each distance
        print(f"For (u,v) = ({u},{v}): r^2 = {u}^2 + ({u})*({v}) + {v}^2 = {r_sq}")
        print(f"-> r = sqrt({r_sq}) = {r:.2f}\n")
        final_distances.append(f"{r:.2f}")

    print("The set of distances is:", final_distances)


if __name__ == "__main__":
    find_planar_distances()