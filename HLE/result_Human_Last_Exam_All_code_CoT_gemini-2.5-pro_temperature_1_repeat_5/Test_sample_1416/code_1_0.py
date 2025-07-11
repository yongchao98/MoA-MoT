import math
import numpy as np

def solve_pyramid_scanning():
    """
    Finds the optimal placement and radii for 6 non-overlapping spherical scans
    inside a pyramid to maximize the total scanned volume.
    """
    # --- Pyramid and Scanner Constants ---
    H = 110.0  # Pyramid height
    L = 75.0   # Pyramid half-base side
    # Pre-calculate constant for plane distance formula
    # Plane equation: 110*|x| + 75*z - 8250 = 0
    # Normalizing factor for the plane vector (110, 0, 75)
    D_SIDE_NORM = math.sqrt(110**2 + 75**2)

    def get_max_radius(x, y, z):
        """
        Calculates the maximum possible radius for a sphere centered at (x, y, z)
        to be fully contained within the pyramid.
        The radius is rounded down to the nearest 0.5m.
        """
        # A point is inside if z is within [0, H] and it's within the side boundaries
        if not (0 <= z <= H and max(abs(x), abs(y)) <= L * (1 - z / H)):
            return 0

        # Constraint 1: Distance to the base plane (z=0)
        max_r_from_base = z
        
        # Constraint 2: Distance to the slanted side planes
        # The point (x,y,z) is inside, so (8250 - 110*max|x|,|y| - 75*z) is positive
        max_r_from_sides = (8250 - 110 * max(abs(x), abs(y)) - 75 * z) / D_SIDE_NORM
        
        # The radius is the minimum of these distances
        max_r = min(max_r_from_base, max_r_from_sides)
        
        # Apply scanner constraints (10m-50m radius, 0.5m steps)
        if max_r > 50.0:
            max_r = 50.0
        if max_r < 10.0:
            return 0
            
        return math.floor(max_r * 2) / 2.0

    # --- Search for Optimal 4+2 Configuration ---
    best_volume_metric = 0
    best_config = None

    # Search space for the middle layer of 4 spheres
    # zm: height of the middle layer
    # d: half the side of the square pattern for the 4 spheres
    # Ranges are chosen based on geometric intuition to be efficient
    for zm in np.arange(20, 35, 0.5):
        for d in np.arange(20, 35, 0.5):
            
            # 1. Calculate properties of the 4 middle spheres
            rm = get_max_radius(d, d, zm)
            if rm == 0: continue
            # Non-overlap condition for adjacent middle spheres: dist(2d) >= r+r
            if d < rm: continue

            # 2. Find optimal bottom sphere on the z-axis
            best_rb = 0
            best_zb = 0
            # Search for zb below the middle layer
            for zb in np.arange(10, zm, 0.5):
                rb_candidate = get_max_radius(0, 0, zb)
                if rb_candidate == 0: continue
                
                # Check for overlap with the middle layer
                dist_sq_bm = 2 * d**2 + (zm - zb)**2
                if dist_sq_bm < (rb_candidate + rm)**2:
                    # Reduce radius if overlapping
                    max_r_by_overlap = math.sqrt(dist_sq_bm) - rm
                    rb_candidate = min(rb_candidate, math.floor(max_r_by_overlap * 2) / 2.0)
                
                if rb_candidate < 10: continue
                
                # We want the largest possible rb for this middle config
                if rb_candidate > best_rb:
                    best_rb = rb_candidate
                    best_zb = zb

            if best_rb == 0: continue
            
            # 3. Find optimal top sphere on the z-axis
            best_rt = 0
            best_zt = 0
            # Search for zt above the middle layer
            for zt in np.arange(zm + rm, H, 0.5):
                rt_candidate = get_max_radius(0, 0, zt)
                if rt_candidate == 0: continue

                # Check overlap with bottom sphere
                if zt - best_zb < rt_candidate + best_rb:
                    max_r_by_overlap_b = (zt - best_zb) - best_rb
                    rt_candidate = min(rt_candidate, math.floor(max_r_by_overlap_b*2)/2.0)
                if rt_candidate < 10: continue

                # Check overlap with middle layer
                dist_sq_tm = 2 * d**2 + (zt - zm)**2
                if dist_sq_tm < (rt_candidate + rm)**2:
                    max_r_by_overlap_m = math.sqrt(dist_sq_tm) - rm
                    rt_candidate = min(rt_candidate, math.floor(max_r_by_overlap_m*2)/2.0)
                if rt_candidate < 10: continue

                if rt_candidate > best_rt:
                    best_rt = rt_candidate
                    best_zt = zt
            
            if best_rt == 0: continue
            
            # 4. Evaluate this complete configuration
            current_volume_metric = 4 * rm**3 + best_rb**3 + best_rt**3
            if current_volume_metric > best_volume_metric:
                best_volume_metric = current_volume_metric
                best_config = {
                    'bottom': {'c': (0, 0, best_zb), 'r': best_rb},
                    'top': {'c': (0, 0, best_zt), 'r': best_rt},
                    'middle_r': rm,
                    'middle_z': zm,
                    'middle_d': d,
                    'radii': [best_rb, best_rt, rm, rm, rm, rm]
                }
    
    # --- Print Results ---
    if best_config:
        print("Optimal scanning configuration found (N=6):")
        print("-" * 40)
        # "Final Equation" here interpreted as the components of the solution
        b = best_config['bottom']
        t = best_config['top']
        rm = best_config['middle_r']
        zm = best_config['middle_z']
        d = best_config['middle_d']
        
        print(f"1. Bottom Sphere: Radius = {b['r']:.1f} m, Center = ({b['c'][0]:.1f}, {b['c'][1]:.1f}, {b['c'][2]:.1f}) m")
        print(f"2. Top Sphere:    Radius = {t['r']:.1f} m, Center = ({t['c'][0]:.1f}, {t['c'][1]:.1f}, {t['c'][2]:.1f}) m")
        print(f"3-6. Middle Layer (4 Spheres):")
        print(f"   Radius = {rm:.1f} m")
        print(f"   Centers at (+/-{d:.1f}, +/-{d:.1f}, {zm:.1f}) m")
        print("-" * 40)

        all_radii = sorted(best_config['radii'])
        r_min = all_radii[0]
        r_max = all_radii[-1]
        
        print(f"Max scanning radius (R): {r_max:.1f} m")
        print(f"Min scanning radius (r): {r_min:.1f} m")
        print("\nAnswer in R:r format:")
        print(f"{r_max:.1f}:{r_min:.1f}")
        
    else:
        print("No valid configuration could be found.")

solve_pyramid_scanning()
<<<25.5:10.0>>>