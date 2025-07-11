import math

def find_planar_distances():
    """
    Calculates the possible normalized distances (r <= 3) between hard spheres
    packed on a 2D plane, considering both hexagonal and square lattice arrangements.
    The distance r is normalized by the sphere diameter d.
    """
    
    # Use a dictionary to store unique distances and their source equation.
    # Key: rounded distance, Value: formatted string showing the calculation.
    # This helps in elegantly handling distances that appear in both lattices (e.g., r=1, 2, 3).
    distance_results = {}

    # Set sphere diameter d=1 for normalized calculations.
    # The search is for r <= 3, which means r^2 <= 9.
    max_r_sq = 9
    
    # A search limit for lattice vector integer components. A limit of 3 is sufficient
    # as 4^2 > 9.
    limit = 3

    print("Deriving distances from possible 2D packing configurations (r <= 3):")
    print("-" * 65)

    # --- Part 1: Hexagonal Lattice ---
    # The squared distance from the origin to a point (n,m) in this lattice
    # is given by r^2 = n^2 + n*m + m^2.
    for n in range(limit + 1):
        for m in range(limit + 1):
            if n == 0 and m == 0:
                continue

            r_sq = float(n**2 + n * m + m**2)
            if r_sq <= max_r_sq + 1e-9: # Use a tolerance for floating point
                r = math.sqrt(r_sq)
                # Store result if this distance is new
                if round(r, 2) not in distance_results:
                    equation = f"r = sqrt({n}^2 + {n}*{m} + {m}^2)"
                    distance_results[round(r, 2)] = f"{equation:<25} = {r:.2f} (from Hexagonal Lattice)"

    # --- Part 2: Square Lattice ---
    # The squared distance from the origin to a point (i,j) in a square lattice
    # is given by r^2 = i^2 + j^2.
    for i in range(limit + 1):
        for j in range(limit + 1):
            # We only need to check one quadrant, e.g., i>=0, j>=0, for unique distances
            if i == 0 and j == 0:
                continue
            
            r_sq = float(i**2 + j**2)
            if r_sq <= max_r_sq + 1e-9: # Use a tolerance
                r = math.sqrt(r_sq)
                # Store result if this distance is new
                if round(r, 2) not in distance_results:
                    equation = f"r = sqrt({i}^2 + {j}^2)"
                    distance_results[round(r, 2)] = f"{equation:<25} = {r:.2f} (from Square Lattice)"

    # --- Part 3: Print the final sorted set of distances ---
    print("\nFinal set of unique distances found:")
    sorted_keys = sorted(distance_results.keys())
    for key in sorted_keys:
        print(distance_results[key])

if __name__ == '__main__':
    find_planar_distances()