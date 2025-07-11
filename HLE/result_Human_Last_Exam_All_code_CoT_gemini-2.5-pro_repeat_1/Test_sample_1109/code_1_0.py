import math
import numpy as np

def solve_seismic_scanning():
    """
    Finds an optimal placement of long-range and short-range seismic scanners
    inside a square pyramid using a greedy, symmetry-aware algorithm.
    """

    # 1. Define constants for the pyramid and scanners
    BASE = 150.0
    HEIGHT = 110.0
    R_LONG = 20.0
    R_SHORT = 7.0
    STEP = 0.5
    
    # Pre-calculate geometric and volume constants for efficiency.
    # The pyramid's side planes are described by equations like 22x + 15z - 1650 = 0.
    # The length of the normal vector to these planes is sqrt(22^2 + 15^2).
    SQRT_709 = math.sqrt(22**2 + 15**2)
    PYRAMID_VOL = (1/3) * (BASE**2) * HEIGHT
    
    # Pre-calculate squared distances for non-overlap checks to avoid costly square roots.
    DIST_SQ_LL = (R_LONG + R_LONG)**2
    DIST_SQ_SS = (R_SHORT + R_SHORT)**2
    DIST_SQ_LS = (R_LONG + R_SHORT)**2

    # 2. Define helper functions for the core logic

    def is_valid_center(center, radius):
        """Checks if a sphere with a given center and radius is fully inside the pyramid."""
        xc, yc, zc = center
        # Condition 1: Sphere must be above the base plane (z=0).
        if zc < radius:
            return False
        # Condition 2: Sphere must be inside the four side planes.
        # This simplifies to a single inequality based on the distance from the center
        # to the nearest side plane.
        if 15 * zc + 22 * max(abs(xc), abs(yc)) > 1650 - radius * SQRT_709:
            return False
        return True

    def get_all_symmetric_points(p):
        """
        Takes a point p=(x,y,z) from the primary octant (x>=y>=0) and returns all
        symmetric points across the axes and y=x plane using a set to handle duplicates.
        """
        px, py, pz = p
        points = set()
        for x_val in {px, -px}:
            for y_val in {py, -py}:
                points.add((x_val, y_val, pz))
                points.add((y_val, x_val, pz))
        return list(points)

    def dist_sq(p1, p2):
        """Calculates the squared Euclidean distance between two points."""
        return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

    # 3. Implement the greedy placement algorithm

    def place_spheres(radius, dist_sq_self, placed_external=None, dist_sq_external=None):
        """A generic function to place spheres of a given radius."""
        placed_spheres = []
        
        # Define the search space for scanner centers in the first octant (x >= y >= 0).
        max_z = (1650 - radius * SQRT_709) / 15
        max_xy_bound = (1650 - radius * SQRT_709 - 15 * radius) / 22
        
        potential_centers = []
        # Use numpy.arange for stable floating-point steps.
        z_coords = np.arange(radius, max_z + STEP, STEP)
        x_coords = np.arange(0, max_xy_bound + STEP, STEP)
        
        for zc in z_coords:
            max_x_at_z = (1650 - radius * SQRT_709 - 15 * zc) / 22
            for xc in x_coords:
                if xc > max_x_at_z: break
                for yc in np.arange(0, xc + STEP, STEP): # Loop y from 0 to x
                    if yc > max_x_at_z: break
                    center = (xc, yc, zc)
                    if is_valid_center(center, radius):
                        potential_centers.append(center)
        
        # Sort potential centers: bottom-up, then center-out for dense packing.
        potential_centers.sort(key=lambda p: (p[2], p[0]**2 + p[1]**2))

        for p in potential_centers:
            sym_points_to_add = get_all_symmetric_points(p)
            can_place = True
            for s_new in sym_points_to_add:
                # Check for overlap with spheres of the same type.
                for s_placed in placed_spheres:
                    if dist_sq(s_new, s_placed) < dist_sq_self:
                        can_place = False; break
                if not can_place: break

                # Check for overlap with spheres of another type (if provided).
                if placed_external:
                    for s_ext in placed_external:
                        if dist_sq(s_new, s_ext) < dist_sq_external:
                            can_place = False; break
                if not can_place: break
            
            if can_place:
                placed_spheres.extend(sym_points_to_add)
        
        return placed_spheres

    # --- Part A: Place Long-Range Scanners ---
    placed_long_spheres = place_spheres(R_LONG, DIST_SQ_LL)

    # --- Part B: Place Short-Range Scanners ---
    placed_short_spheres = place_spheres(R_SHORT, DIST_SQ_SS, 
                                         placed_external=placed_long_spheres, 
                                         dist_sq_external=DIST_SQ_LS)

    # 4. Calculate final metrics and print the result.
    n = len(placed_long_spheres)
    m = len(placed_short_spheres)
    
    vol_long = n * (4/3) * math.pi * (R_LONG**3)
    vol_short = m * (4/3) * math.pi * (R_SHORT**3)
    total_scanned_volume = vol_long + vol_short
    
    coverage_ratio = (total_scanned_volume / PYRAMID_VOL) * 100
    
    # Output the final numbers in the specified format.
    print(f"{n}:{m}:{coverage_ratio:.1f}%")

# Execute the main function
solve_seismic_scanning()
<<<13:84:48.7%>>>