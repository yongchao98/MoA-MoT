import math

def find_hard_sphere_distances():
    """
    Calculates and prints the possible normalized distances (r) between
    hard spheres on a 2D plane for r <= 3.

    The distances are derived from the geometry of square and hexagonal lattices,
    which represent the fundamental packing arrangements.
    """
    r_max = 3.0
    # Use a dictionary to store unique distances and their calculation string
    # Keys are the distances rounded to handle floating point inaccuracies
    # Values are tuples of (exact_distance, calculation_string)
    found_distances = {}

    # Iterate through possible integer coordinates (n, m)
    # A search limit of int(r_max) is sufficient for n and m.
    limit = int(r_max)
    for n in range(limit + 1):
        for m in range(limit + 1):
            if n == 0 and m == 0:
                continue

            # --- Case 1: Square Lattice ---
            # r^2 = n^2 + m^2
            r_sq_square = float(n**2 + m**2)
            if r_sq_square <= r_max**2:
                r_square = math.sqrt(r_sq_square)
                # Use a rounded key to avoid floating point duplicates
                key = round(r_square, 5)
                if key not in found_distances:
                    # To make the equation simpler, we show the non-zero terms
                    if n == 0 or m == 0:
                        term = max(n, m)
                        calc_string = f"r = sqrt({term}^2)"
                    else:
                        calc_string = f"r = sqrt({n}^2 + {m}^2)"
                    found_distances[key] = (r_square, calc_string)

            # --- Case 2: Hexagonal Lattice ---
            # r^2 = n^2 + m^2 + nm
            r_sq_hex = float(n**2 + m**2 + n * m)
            if r_sq_hex <= r_max**2:
                r_hex = math.sqrt(r_sq_hex)
                key = round(r_hex, 5)
                if key not in found_distances:
                    if n == 0 or m == 0:
                        term = max(n, m)
                        calc_string = f"r = sqrt({term}^2)"
                    else:
                        calc_string = f"r = sqrt({n}^2 + {m}^2 + {n}*{m})"
                    found_distances[key] = (r_hex, calc_string)

    # Sort the results based on the distance
    sorted_distances = sorted(found_distances.values(), key=lambda item: item[0])

    print("The set of possible distances r for r <= 3 are:")
    for r, calc in sorted_distances:
        print(f"{calc} = {r:.2f}")

find_hard_sphere_distances()