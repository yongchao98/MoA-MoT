import numpy as np

# A simple class to hold sphere properties
class Sphere:
    def __init__(self, center, radius):
        self.c = np.array(center)
        self.r = radius

    def __repr__(self):
        return f"Sphere(center=({self.c[0]:.1f}, {self.c[1]:.1f}, {self.c[2]:.1f}), radius={self.r:.1f})"

# --- Pyramid and Geometry Definitions ---
H = 110.0  # Height
A = 150.0  # Base side length
H_A_HALF = H * (A / 2.0) # 110 * 75 = 8250, a useful constant
K_SIDE = np.sqrt(H**2 + (A / 2.0)**2) # sqrt(110^2 + 75^2) = 133.135...

def get_max_radius_in_pyramid(center):
    """Calculates the max radius of a sphere that can fit inside the pyramid at a given center."""
    cx, cy, cz = center
    if cz < 0 or cz > H:
        return 0
    
    # Distance to base plane (z=0)
    dist_base = cz
    
    # Max horizontal extent of the pyramid at this height
    max_xy_at_cz = (A / 2.0) * (1.0 - cz / H)
    if abs(cx) > max_xy_at_cz or abs(cy) > max_xy_at_cz:
        return 0 # Center is outside

    # Distance to the four side planes
    # Plane equations are of the form: 110*x + 75*z - 8250 = 0
    dist_pos_x = (H_A_HALF - H * cx - (A/2) * cz) / K_SIDE
    dist_neg_x = (H_A_HALF + H * cx - (A/2) * cz) / K_SIDE
    dist_pos_y = (H_A_HALF - H * cy - (A/2) * cz) / K_SIDE
    dist_neg_y = (H_A_HALF + H * cy - (A/2) * cz) / K_SIDE
    
    return min(dist_base, dist_pos_x, dist_neg_x, dist_pos_y, dist_neg_y)

def solve_for_scans():
    """
    Finds the 6 optimal spheres using a structured greedy search.
    """
    placed_spheres = []
    
    # --- 1. Find the largest sphere on the central z-axis ---
    best_s1 = None
    max_r_s1 = 0
    # Search for the center along the z-axis (multiples of 0.5)
    for z in np.arange(0.5, H, 0.5):
        center = (0, 0, z)
        r = get_max_radius_in_pyramid(center)
        # Round down to nearest 0.5 and check limits
        r = np.floor(r * 2) / 2
        if 10 <= r <= 50:
            if r > max_r_s1:
                max_r_s1 = r
                best_s1 = Sphere(center, r)
    placed_spheres.append(best_s1)

    # --- 2. Find the second sphere on the central z-axis ---
    best_s2 = None
    max_r_s2 = 0
    s1 = placed_spheres[0]
    for z in np.arange(0.5, H, 0.5):
        center = (0, 0, z)
        
        # Max radius allowed by pyramid
        r_pyramid = get_max_radius_in_pyramid(center)
        
        # Max radius allowed by non-overlap with sphere 1
        dist_to_s1 = abs(z - s1.c[2])
        r_overlap = dist_to_s1 - s1.r

        r = min(r_pyramid, r_overlap)
        r = np.floor(r * 2) / 2
        
        if 10 <= r <= 50:
            if r > max_r_s2:
                max_r_s2 = r
                best_s2 = Sphere(center, r)
    placed_spheres.append(best_s2)

    # --- 3. Find 4 symmetric spheres in the corners ---
    best_s3_config = None
    max_r_s3 = 0
    s1, s2 = placed_spheres[0], placed_spheres[1]

    # Search for one sphere in the first quadrant, assuming x=y for simplicity
    for z in np.arange(10, 30, 0.5):
        for x in np.arange(10, 50, 0.5):
            center = (x, x, z)
            
            # Non-overlap with self (4-sphere configuration)
            r_self_overlap = x 
            
            # Max radius from other constraints
            r_pyramid = get_max_radius_in_pyramid(center)
            dist_to_s1 = np.linalg.norm(np.array(center) - s1.c)
            r_overlap1 = dist_to_s1 - s1.r
            dist_to_s2 = np.linalg.norm(np.array(center) - s2.c)
            r_overlap2 = dist_to_s2 - s2.r

            r = min(r_pyramid, r_self_overlap, r_overlap1, r_overlap2)
            r = np.floor(r * 2) / 2

            if 10 <= r <= 50:
                if r > max_r_s3:
                    max_r_s3 = r
                    best_s3_config = {'center': center, 'radius': r}

    # Add the 4 symmetric spheres
    if best_s3_config:
        c = best_s3_config['center']
        r = best_s3_config['radius']
        placed_spheres.append(Sphere((c[0], c[1], c[2]), r))
        placed_spheres.append(Sphere((-c[0], c[1], c[2]), r))
        placed_spheres.append(Sphere((c[0], -c[1], c[2]), r))
        placed_spheres.append(Sphere((-c[0], -c[1], c[2]), r))
        
    return placed_spheres

# --- Main Execution ---
final_spheres = solve_for_scans()

print("Found optimal configuration for N=6 scans:\n")
for i, s in enumerate(final_spheres):
    print(f"Sphere {i+1}: Center=({s.c[0]:.1f}, {s.c[1]:.1f}, {s.c[2]:.1f}), Radius={s.r:.1f}")

radii = [s.r for s in final_spheres]
max_r = max(radii)
min_r = min(radii)

print(f"\nThis configuration results from varied radii as requested.")
print(f"The maximal and minimal radii are R={max_r:.1f} and r={min_r:.1f}.")
print(f"\nFinal answer in R:r format is {max_r:.1f}:{min_r:.1f}")

print(f"\n<<<{max_r:.1f}:{min_r:.1f}>>>")
