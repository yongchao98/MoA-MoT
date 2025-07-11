import math

def find_hard_sphere_distances():
    """
    Calculates the possible normalized distances (r) between hard spheres
    in a 2D jammed packing for r <= 3.

    The model assumes a triangular lattice, which is the densest packing in 2D.
    The squared distance r^2 is given by the formula n^2 + nm + m^2, where
    n and m are integer multiples of the lattice basis vectors.
    """

    # Use a set to store unique rounded distances to avoid printing duplicates.
    found_distances_rounded = set()
    # Use a list to store results as (distance, n, m) for sorted output.
    results = []
    
    # Iterate through a reasonable range of integer multipliers n and m.
    # If |n| or |m| > 3, r will likely be > 3. We use a limit of 4 to be safe.
    limit = 4

    # We can iterate with n>=0 due to the symmetry of the lattice.
    for n in range(0, limit + 1):
        for m in range(-limit, limit + 1):
            
            # Calculate the squared normalized distance
            r_squared = n**2 + n * m + m**2

            # Skip the origin point (distance to self)
            if r_squared == 0:
                continue

            # Calculate the distance
            r = math.sqrt(r_squared)

            # Check if the distance is within the desired range (r <= 3)
            # Use a small epsilon for floating-point comparison
            if r <= 3.00001:
                # Format the distance to two decimal places for comparison
                r_rounded_str = f"{r:.2f}"
                
                # If this distance is new, store it
                if r_rounded_str not in found_distances_rounded:
                    found_distances_rounded.add(r_rounded_str)
                    results.append((r, n, m))

    # Sort the results by distance
    results.sort()

    print("The set of possible normalized distances r is found using the equation r = sqrt(n^2 + n*m + m^2):")
    print("-" * 60)
    for r_val, n_val, m_val in results:
        # Fulfills the requirement to "output each number in the final equation"
        print(f"For (n,m) = ({n_val},{m_val}): r = sqrt({n_val}^2 + {n_val}*{m_val} + {m_val}^2) = {r_val:.2f}")
    
    print("-" * 60)
    print("The final unique distances are:")
    # Print the sorted, unique values
    final_list_sorted = sorted([float(d) for d in found_distances_rounded])
    for dist in final_list_sorted:
        print(f"{dist:.2f}")

find_hard_sphere_distances()
