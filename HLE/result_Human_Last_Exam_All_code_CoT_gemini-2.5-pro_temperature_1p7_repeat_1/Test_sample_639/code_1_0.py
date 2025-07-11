import math

def find_planar_distances():
    """
    Calculates the possible distances between centers of hard spheres
    packed on a 2D plane for normalized distances r <= 3.
    """
    
    # Use a dictionary to store the results, keyed by the rounded distance
    # to handle floating point inaccuracies and ensure uniqueness.
    # The value will be a tuple: (exact_r, n, m, r_squared_int)
    results = {}

    # Iterate through a reasonable range of integer coefficients n and m.
    # Iterating with n>=0 and a full range for m is sufficient to find
    # all unique distances due to the symmetries of the lattice.
    for n in range(0, 4):
        for m in range(-4, 4):
            # Skip the case (0,0) which is the distance to the particle itself
            if n == 0 and m == 0:
                continue

            # Calculate the squared distance using the derived formula
            r_squared = n**2 + m**2 + n*m
            
            # Avoid taking sqrt of negative number if by some chance r_squared is negative
            if r_squared < 0:
                continue

            r = math.sqrt(r_squared)
            
            # Check if the distance is within the specified range 0 < r <= 3.
            # Add a small epsilon for floating point comparison.
            if 1e-9 < r <= 3.00001:
                # Use the two-decimal rounded value as a key to group results
                r_rounded_key = round(r, 2)
                
                # If we haven't found this distance yet, store it with its calculation details.
                if r_rounded_key not in results:
                    # Store the exact distance r, the (n, m) pair, and integer r_squared
                    results[r_rounded_key] = (r, n, m, r_squared)

    # Get the keys (which are the rounded distances) and sort them
    sorted_keys = sorted(results.keys())

    print("The set of possible normalized distances r <= 3 for hard spheres on a plane are found using the equation r = sqrt(n^2 + m^2 + n*m).")
    print("The distinct distances and their origins are:")
    
    # Print the results in a formatted way
    for key in sorted_keys:
        r_exact, n, m, r_sq_int = results[key]
        # Use a proper Unicode symbol for approx equals
        print(f"r = sqrt({n}\u00b2 + {m}\u00b2 + {n}*{m}) = sqrt({r_sq_int}) \u2248 {r_exact:.2f}")

# Execute the function to find and print the distances
find_planar_distances()