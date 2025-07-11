import numpy as np

def solve_pyramid_scanning():
    """
    Finds the optimal radii for 6 scans in a 3x2 grid configuration inside a pyramid.
    """
    # Pyramid dimensions
    base_side = 150.0
    height = 110.0
    half_base = base_side / 2

    # Scanner constraints
    min_radius = 10.0
    max_radius = 50.0
    step = 0.5

    # Brute-force search for the optimal radii combination
    best_volume = 0
    best_config = None
    
    radii_range = np.arange(min_radius, max_radius + step, step)

    # Iterate through all possible radii for outer (ro) and central (rc) spheres
    for ro in radii_range:
        for rc in radii_range:
            # Determine grid parameters based on tight packing
            zc = max(ro, rc)
            y0 = max(ro, rc)
            d = ro + rc

            # Constraint for the corner sphere (radius ro)
            # Center is at (d, y0, zc). Check if this sphere is inside.
            # The condition is max(|xc|,|yc|) <= (a/2)*(1-(zc+r)/h)
            # max(|d|,|y0|) = max(ro+rc, max(ro,rc)) = ro+rc
            corner_constraint_ok = (ro + rc) <= half_base * (1 - (zc + ro) / height)

            # Constraint for the central sphere (radius rc)
            # Center is at (0, y0, zc). Check if this sphere is inside.
            # max(|xc|,|yc|) = y0 = max(ro, rc)
            center_constraint_ok = y0 <= half_base * (1 - (zc + rc) / height)
            
            if corner_constraint_ok and center_constraint_ok:
                # If both spheres fit, calculate the total volume proxy
                current_volume = 4 * ro**3 + 2 * rc**3
                if current_volume > best_volume:
                    best_volume = current_volume
                    best_config = {
                        'ro': ro, 'rc': rc, 'zc': zc, 
                        'd': d, 'y0': y0, 'volume': current_volume
                    }

    # Output the results
    if best_config:
        ro = best_config['ro']
        rc = best_config['rc']
        zc = best_config['zc']
        d = best_config['d']
        y0 = best_config['y0']
        
        print("Optimal Configuration Found:")
        print("-" * 30)
        print(f"Configuration type: 3x2 planar grid")
        print(f"Grid is placed at height zc = {zc:.1f} m")
        print("\n--- Sphere Details ---")
        print(f"4 Corner Spheres:")
        print(f"  Radius (ro): {ro:.1f} m")
        print(f"  Centers (x,y,z): (±{d:.1f}, ±{y0:.1f}, {zc:.1f})")
        print(f"\n2 Central Spheres:")
        print(f"  Radius (rc): {rc:.1f} m")
        print(f"  Centers (x,y,z): (0, ±{y0:.1f}, {zc:.1f})")
        print("-" * 30)

        max_r = max(ro, rc)
        min_r = min(ro, rc)
        print(f"\nMax scanning radius (R): {max_r:.1f} m")
        print(f"Min scanning radius (r): {min_r:.1f} m")
        
        # Final answer format
        print("\nFinal Answer:")
        print(f"{max_r:.1f}:{min_r:.1f}")
        
        return f"<<<{max_r:.1f}:{min_r:.1f}>>>"

    else:
        print("No valid configuration was found.")
        return "<<<No solution found>>>"

# Execute the function and capture the final answer string
final_answer_str = solve_pyramid_scanning()
# The final answer is wrapped for the system to extract.
# print(final_answer_str)