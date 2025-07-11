import math

def calculate_area():
    """
    Calculates the area of the region R defined by floor(|z|) = |floor(z)|
    within the square [0,6]x[0,6] using numerical integration.
    """
    grid_per_unit = 2000  # Grid resolution for numerical integration
    step = 1.0 / grid_per_unit
    delta_area = step * step

    # Step 1 & 2: Identify valid unit squares (n,m) where k = sqrt(n^2+m^2) is an integer.
    valid_pairs = []
    for n in range(6):
        for m in range(6):
            val = n**2 + m**2
            k = math.isqrt(val)
            if k**2 == val:
                valid_pairs.append(((n, m), k))

    # Step 3 & 4: Calculate area for each unique square shape and store in a cache.
    # Due to symmetry, Area(n,m) = Area(m,n). We cache results for sorted pairs (min(n,m), max(n,m)).
    area_cache = {}
    
    # A dictionary to hold the components of the final sum for printing.
    equation_components = {}

    for (n, m), k in valid_pairs:
        # Use a sorted tuple as a key to handle symmetry, e.g. (3,4) and (4,3) are the same shape.
        pair_key = tuple(sorted((n, m)))
        if pair_key not in area_cache:
            # Calculate the area for this square shape if not already computed.
            square_area = 0
            k_sq = k**2
            k_plus_1_sq = (k + 1)**2
            # For the numerical integration, we iterate over a fine grid within the unit square.
            # We use the square at (n,m) for the calculation.
            for i in range(grid_per_unit):
                x = n + (i + 0.5) * step
                for j in range(grid_per_unit):
                    y = m + (j + 0.5) * step
                    
                    z_mod_sq = x**2 + y**2
                    
                    # Check if the point is in the required annulus: k^2 <= |z|^2 < (k+1)^2
                    if k_sq <= z_mod_sq < k_plus_1_sq:
                        square_area += delta_area
            
            area_cache[pair_key] = square_area

    # Step 5: Sum up the areas from all valid squares to get the total area.
    total_area = 0.0
    for (n, m), k in valid_pairs:
        pair_key = tuple(sorted((n, m)))
        total_area += area_cache[pair_key]
    
    print("The total area is a sum of contributions from unit squares [n,n+1)x[m,m+1) where sqrt(n^2+m^2) is an integer.")
    print("The individual area contributions are calculated numerically.")
    print("\nEquation for the total area:")

    # Prepare strings for the final equation output
    sum_parts = []
    for pair_key in sorted(area_cache.keys()):
      area = area_cache[pair_key]
      n, m = pair_key
      if n == m:
        # Diagonal squares like (0,0) appear once
        sum_parts.append(f"{area:.4f}")
      else:
        # Off-diagonal squares like (0,1) and (1,0) appear twice by symmetry
        sum_parts.append(f"2*{area:.4f}")

    print(f"Total Area = {' + '.join(sum_parts)}")
    print(f"\nResult:")
    print(f"The total area is {total_area:.2f}")

calculate_area()