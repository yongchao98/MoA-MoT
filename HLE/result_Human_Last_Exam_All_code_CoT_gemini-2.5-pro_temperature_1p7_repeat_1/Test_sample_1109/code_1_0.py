import math
import numpy as np
import sys

def solve_seismic_scanning():
    """
    Solves the seismic scanning puzzle by finding an optimal placement
    of long-range and short-range scanners inside a pyramid.
    """
    # 1. Define constants and helper functions
    B = 150.0  # Pyramid base side
    H = 110.0  # Pyramid height
    r_L = 20.0 # Long range radius
    r_S = 7.0  # Short range radius
    grid_step = 0.5

    # Volume calculations
    V_pyramid = (1/3) * B**2 * H
    V_sphere_L = (4/3) * math.pi * r_L**3
    V_sphere_S = (4/3) * math.pi * r_S**3

    def is_valid_center(x, y, z, r, B, H):
        """
        Checks if a sphere of radius r centered at (x,y,z) is
        entirely within the pyramid.
        The center must be within an 'inner' pyramid whose faces are
        shrunk inwards by a distance r.
        """
        # Condition 1: Sphere must be above the base plane (z=0)
        if z < r:
            return False

        # Condition 2: Sphere must be inside the four slanted side faces
        # This is checked by ensuring the center point is far enough from each face.
        # The distance check is simplified by projecting it onto the xy plane.
        # A sphere is inside if its center is within an inner pyramid.
        d_slant = r * math.sqrt(4 * H**2 + B**2)
        if (2 * H * abs(x) + B * z > B * H - d_slant) or \
           (2 * H * abs(y) + B * z > B * H - d_slant):
            return False
        
        return True

    # 2. Greedy placement algorithm
    scanners = []
    radii_to_place = [r_L, r_S]

    for r in radii_to_place:
        print(f"\n[+] Placing scanners with radius {r}m...", file=sys.stderr)
        
        # Define search space bounds for this radius
        d_slant = r * math.sqrt(4 * H**2 + B**2)
        z_max_center = H - (d_slant / B)
        
        z_coords = np.arange(r, z_max_center + grid_step, grid_step)

        # Iterate z from bottom to top
        for z in z_coords:
            # Determine valid x/y range for this height z
            max_coord_at_z = (B*H - d_slant - B*z) / (2*H)
            
            # Create a 2D grid of candidate points for the current z-slice
            # and sort them to start from the center and expand outwards.
            # This is a key part of the greedy strategy.
            x_range = np.arange(0, max_coord_at_z + grid_step, grid_step)
            y_range = np.arange(0, max_coord_at_z + grid_step, grid_step)
            
            coords_2d = []
            for x in x_range:
                for y in y_range:
                    if x == 0 and y == 0:
                        coords_2d.append((0, 0))
                    elif x == 0:
                        coords_2d.append((0, y)); coords_2d.append((0, -y))
                    elif y == 0:
                        coords_2d.append((x, 0)); coords_2d.append((-x, 0))
                    else:
                        coords_2d.append((x, y)); coords_2d.append((-x, y))
                        coords_2d.append((x, -y)); coords_2d.append((-x, -y))
            
            # Use set to get unique coordinates and sort by distance from origin
            unique_coords = sorted(list(set(coords_2d)), key=lambda p: p[0]**2 + p[1]**2)

            for x, y in unique_coords:
                candidate_center = (x, y, z)
                
                # We know the center is valid due to the loop bounds.
                # Now check for overlap with already placed scanners.
                is_overlapped = False
                for s in scanners:
                    s_center = s['center']
                    dist_sq = (candidate_center[0] - s_center[0])**2 + \
                              (candidate_center[1] - s_center[1])**2 + \
                              (candidate_center[2] - s_center[2])**2
                    min_dist_sq = (r + s['radius'])**2
                    
                    if dist_sq < min_dist_sq - 1e-9: # a small tolerance for float precision
                        is_overlapped = True
                        break
                
                if not is_overlapped:
                    scanners.append({'center': candidate_center, 'radius': r})

    print(f"\n[+] Placement complete.", file=sys.stderr)

    # 3. Calculate and report results
    n = sum(1 for s in scanners if s['radius'] == r_L)
    m = sum(1 for s in scanners if s['radius'] == r_S)

    vol_from_L = n * V_sphere_L
    vol_from_S = m * V_sphere_S
    total_scanned_volume = vol_from_L + vol_from_S
    coverage = total_scanned_volume / V_pyramid

    print("\n--- Scan Analysis ---")
    print(f"Pyramid Volume: {V_pyramid:,.2f} m^3")
    print("-" * 25)
    print(f"Long range scans (n): {n}")
    print(f"  - Volume from long range: {n} * {V_sphere_L:,.2f} = {vol_from_L:,.2f} m^3")
    print(f"Short range scans (m): {m}")
    print(f"  - Volume from short range: {m} * {V_sphere_S:,.2f} = {vol_from_S:,.2f} m^3")
    print("-" * 25)
    print(f"Total Scanned Volume: {total_scanned_volume:,.2f} m^3")
    print(f"Coverage Ratio (p): {coverage:.3%} (rounded to {coverage:.1%})")
    print("\n--- Final Answer ---")
    
    # Final answer format as requested n:m:p
    final_answer = f"{n}:{m}:{coverage*100:.1f}"
    print(final_answer)
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    solve_seismic_scanning()